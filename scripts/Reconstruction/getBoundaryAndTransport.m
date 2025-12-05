function [exchangeRxns, exchangeRxnsIndexes, pureTransportRxns, pureTransportIndexes] = getBoundaryAndTransport(model, reactionType)
% getBoundaryAndTransport
%   Unified detector for:
%     (A) Exchange/boundary reactions (same API options as getExchangeRxns)
%     (B) “Pure same-entity compartment switch” transport reactions
%         (exactly two metabolites, opposite & equal stoich, different compartments,
%          and identical Name/Formula/Charge).
%
% Usage:
%   [exRxns, exIdx, ptRxns, ptIdx] = getBoundaryAndTransport(model, 'all')
%   [~, ~, ptRxns, ptIdx]          = getBoundaryAndTransport(model, 'pureTransport')   % or 'pure'
%
% Inputs:
%   model        : COBRA-like model (requires S, rxns, mets; optionally metNames, metFormulas,
%                  metCharges, metComps, unconstrained, lb, ub)
%   reactionType : For EXCHANGE/BORDER subset (first two outputs), choose one of:
%                  'all' | 'in' | 'out' | 'blocked' | 'reverse' | 'uptake' | 'excrete'
%                  Additionally supports: 'pureTransport' (alias 'pure')
%                  - If 'pureTransport'/'pure' is passed, the exchange outputs are returned empty
%                    and pure-transport outputs are still returned as usual.
%
% Outputs:
%   exchangeRxns         : cellstr of exchange/boundary reaction IDs (filtered by reactionType)
%   exchangeRxnsIndexes  : column vector of indices of those reactions
%   pureTransportRxns    : cellstr of “pure same-entity compartment switch” reaction IDs (always computed)
%   pureTransportIndexes : column vector of indices of those reactions
%
% Notes:
%   - The function ALWAYS computes pureTransport detection and returns it (3rd, 4th outputs),
%     regardless of reactionType.
%   - You can call with 'all' in Step 3 to exclude both exchange and pureTransport from irreversible splitting:
%       [~, exIdx, ~, ptIdx] = getBoundaryAndTransport(model, 'all');
%       toSplitMask(exIdx) = false; toSplitMask(ptIdx) = false;
%   - Calling with 'pureTransport' is useful when you explicitly want to operate only on that set.

    if nargin < 2 || isempty(reactionType)
        reactionType = 'all';
    else
        reactionType = char(reactionType);
    end

    nRxn = numel(model.rxns);

    %% ---------- (A) Exchange/boundary detection (compatible with getExchangeRxns) ----------
    hasNoProd = false(nRxn,1);
    hasNoSubs = false(nRxn,1);

    if isfield(model, 'unconstrained') && ~isempty(model.unconstrained)
        rowSel = (model.unconstrained ~= 0);
        if any(rowSel)
            [~, Ipos] = find(model.S(rowSel, :) > 0);
            hasNoProd(Ipos) = true;
            [~, Ineg] = find(model.S(rowSel, :) < 0);
            hasNoSubs(Ineg) = true;
        end
    else
        % Column-wise logic (per reaction): no products if no positive stoich, no substrates if no negative stoich
        hasNoProd = (sum(model.S > 0, 1) == 0).';
        hasNoSubs = (sum(model.S < 0, 1) == 0).';
    end

    allExchIdx = find(hasNoProd | hasNoSubs);  % indices for 'all' exchange/boundary

    switch lower(reactionType)
        case {'both','all'}
            exchangeRxnsIndexes = allExchIdx;

        case 'in'
            exchangeRxnsIndexes = find(hasNoSubs);

        case 'out'
            exchangeRxnsIndexes = find(hasNoProd);

        case 'blocked'
            mask = false(nRxn,1);
            mask(allExchIdx) = (model.lb(allExchIdx) == 0 & model.ub(allExchIdx) == 0);
            exchangeRxnsIndexes = find(mask);

        case 'reverse'
            mask = false(nRxn,1);
            mask(allExchIdx) = (model.lb(allExchIdx) < 0 & model.ub(allExchIdx) > 0);
            exchangeRxnsIndexes = find(mask);

        case 'uptake'
            mask = false(nRxn,1);
            up1 = find(hasNoProd & (model.lb < 0) & (model.ub <= 0));
            up2 = find(hasNoSubs & (model.lb >= 0) & (model.ub > 0));
            mask([up1; up2]) = true;
            exchangeRxnsIndexes = find(mask);

        case 'excrete'
            mask = false(nRxn,1);
            ex1 = find(hasNoProd & (model.lb >= 0) & (model.ub > 0));
            ex2 = find(hasNoSubs & (model.lb < 0)  & (model.ub <= 0));
            mask([ex1; ex2]) = true;
            exchangeRxnsIndexes = find(mask);

        case {'puretransport','pure'}
            % Caller explicitly wants pure-transport set; make exchange outputs empty for clarity.
            exchangeRxnsIndexes = zeros(0,1);

        otherwise
            error('Invalid reactionType specified: %s', reactionType);
    end

    exchangeRxnsIndexes = sort(exchangeRxnsIndexes(:));
    exchangeRxns        = model.rxns(exchangeRxnsIndexes);

    %% ---------- (B) Pure same-entity compartment switch detection ----------
    pureTransportMask = false(nRxn,1);

    hasNames   = isfield(model,'metNames')    && numel(model.metNames)    == numel(model.mets);
    hasForm    = isfield(model,'metFormulas') && numel(model.metFormulas) == numel(model.mets);
    hasCharge  = isfield(model,'metCharges')  && numel(model.metCharges)  == numel(model.mets);
    hasMetComp = isfield(model,'metComps')    && numel(model.metComps)    == numel(model.mets);

    for j = 1:nRxn
        nz = find(model.S(:,j) ~= 0);
        if numel(nz) ~= 2, continue; end

        s = full(model.S(nz, j));
        if ~(s(1) == -s(2)), continue; end   % opposite & equal magnitude

        % compartments must differ
        if hasMetComp
            c1 = model.metComps(nz(1)); c2 = model.metComps(nz(2));
            if c1 == c2, continue; end
        else
            comp1 = parseCompFromMet(model.mets{nz(1)});
            comp2 = parseCompFromMet(model.mets{nz(2)});
            if ~isempty(comp1) && strcmp(comp1, comp2), continue; end
        end

        % chemical identity equality (when each field is available)
        sameName   = true;
        sameForm   = true;
        sameCharge = true;

        if hasNames
            sameName = strcmp(model.metNames{nz(1)}, model.metNames{nz(2)});
        end
        if hasForm
            f1 = model.metFormulas{nz(1)}; f2 = model.metFormulas{nz(2)};
            sameForm = ~(isempty(f1) ~= isempty(f2)) && strcmp(f1, f2);
        end
        if hasCharge
            sameCharge = (model.metCharges(nz(1)) == model.metCharges(nz(2)));
        end

        if sameName && sameForm && sameCharge
            pureTransportMask(j) = true;
        end
    end

    pureTransportIndexes = find(pureTransportMask);
    pureTransportRxns    = model.rxns(pureTransportIndexes);
end

% --- helper: parse compartment suffix like "glc__D[c]" into "c"
function comp = parseCompFromMet(metId)
    tok = regexp(metId, '\[([^\]]+)\]$', 'tokens','once');
    if isempty(tok), comp = ''; else, comp = tok{1}; end
end
