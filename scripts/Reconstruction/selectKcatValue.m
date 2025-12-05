function [model, rxnIdx] = selectKcatValue(model, kcatList, criteria, overwrite, opts)
% selectKcatValue (ECOMAP, prioritized)
%   Select ONE kcat per reaction from a pooled candidate list using
%   (1) source-aware priority, (2) fuzzy specificity (origin/wildcardLvl),
%   and (3) a statistical tie-breaker (max|min|median|mean).
%   Writes the chosen value into model.enzymeConstraints.kcat and tracks
%   provenance in model.enzymeConstraints.source.
%
%   ─────────────────────────────────────────────────────────────────────
%   Origin & Attribution
%     • This routine is adapted from the GECKO 3 codebase (selectKcatValue),
%       then extended for ECOMAP to leverage CompleteMatch vs. Fuzzy
%       (origin/wildcards) vs. Prediction sources with configurable
%       priority and robust tie-breaking.
%   ─────────────────────────────────────────────────────────────────────
%
% Inputs
%   model     : ECOMAP model; requires fields
%                 - model.enzymeConstraints.rxns  (cellstr, N×1)
%                 - model.enzymeConstraints.kcat  (double, N×1)      [created if missing]
%                 - model.enzymeConstraints.source(cellstr, N×1)     [created if missing]
%   kcatList  : struct of candidate entries (one row per candidate) with:
%                 - rxns        : reaction IDs (must exist in enzymeConstraints.rxns)
%                 - kcats       : numeric kcat (1/s)
%                 - kcatSource  : per-row provenance (e.g., 'CompleteMatch','brenda','DLKcat',...)
%                                  (If absent but .source exists, it is expanded into per-row kcatSource)
%                 - genes       : (optional) gene IDs; rows having genes are preferred if ties remain
%                 - substrates  : (optional) substrate names
%                 - origin      : (optional, for fuzzy) specificity tier (1..6, smaller is better)
%                 - wildcardLvl : (optional, for fuzzy) EC wildcard level (0..3, smaller is better)
%   criteria  : tie-break rule within the final candidate tier for each reaction:
%                 'max' (default) | 'min' | 'median' | 'mean'
%   overwrite : write policy for existing kcat values:
%                 'true' (default)  → always overwrite
%                 'false'           → only write where current kcat == 0
%                 'ifHigher'        → overwrite only if new > current
%               (logical true/false also accepted)
%   opts      : optional struct to customize behavior
%                 .priority          : cellstr high→low source order
%                       default {'CompleteMatch','brenda','DLKcat','UniKP','CatPred'}
%                 .topOriginLimit    : exact-EC fuzzy with origin ≤ this are high-priority (default 6)
%                 .bottomOriginLimit : lower-priority fuzzy allowed up to this origin (default 6)
%                 .wildcardLimit     : lower-priority fuzzy allowed up to this wildcard level (default 3)
%                 .capKcat           : numeric upper cap for kcat (diffusion limit; default 1e7)
%
% Outputs
%   model  : updated model with chosen kcat/source written at affected indices
%   rxnIdx : indices (into model.enzymeConstraints.rxns) that were updated
%
% Example
%   [model, idx] = selectKcatValue(model, kcatList_merged, 'max', 'ifHigher');

    % ------------------------ defaults ----------------------------------
    if nargin < 5 || isempty(opts), opts = struct(); end
    if nargin < 4 || isempty(overwrite), overwrite = 'true'; end
    if islogical(overwrite), overwrite = tern(overwrite,'true','false'); end
    if nargin < 3 || isempty(criteria), criteria = 'max'; end

    if ~isfield(opts,'priority') || isempty(opts.priority)
        opts.priority = {'CompleteMatch','brenda','DLKcat','UniKP','CatPred'};
    end
    if ~isfield(opts,'topOriginLimit'),     opts.topOriginLimit = 6; end
    if ~isfield(opts,'bottomOriginLimit'),  opts.bottomOriginLimit = 6; end
    if ~isfield(opts,'wildcardLimit'),      opts.wildcardLimit = 3; end
    if ~isfield(opts,'capKcat'),            opts.capKcat = 1e7; end

    ec = model.enzymeConstraints;
    nRxn = numel(ec.rxns);
    if ~isfield(ec,'kcat')   || numel(ec.kcat)   ~= nRxn, ec.kcat   = zeros(nRxn,1);   end
    if ~isfield(ec,'source') || numel(ec.source) ~= nRxn, ec.source = cell(nRxn,1);    end

    % Ensure per-row kcatSource is available
    if ~isfield(kcatList,'kcatSource') || isempty(kcatList.kcatSource)
        if isfield(kcatList,'source')
            kcatList.kcatSource = repmat({kcatList.source}, numel(kcatList.kcats), 1);
        else
            error('kcatList must contain either .kcatSource (per-row) or .source (scalar).');
        end
    end

    % ------------------- basic cleaning & capping -----------------------
    k = kcatList.kcats;
    k(~isfinite(k))  = NaN;           % drop non-finite
    k(k <= 0)        = NaN;           % drop non-positive
    k(k > opts.capKcat) = opts.capKcat; % cap to diffusion limit
    bad = isnan(k);

    % Remove bad rows across aligned fields
    fields = {'kcats','rxns','kcatSource','genes','substrates','origin','wildcardLvl','eccodes'};
    for f = 1:numel(fields)
        if isfield(kcatList, fields{f})
            tmp = kcatList.(fields{f});
            if numel(tmp) == numel(bad)
                kcatList.(fields{f})(bad) = [];
            end
        end
    end
    if isempty(kcatList.kcats)
        rxnIdx = [];
        model.enzymeConstraints = ec;
        return;
    end

    % -------------------- map to model rxn indices ----------------------
    [ok, idxInModel] = ismember(kcatList.rxns, ec.rxns);
    if ~all(ok)
        warning('Some kcatList.rxns are not present in model.enzymeConstraints.rxns.');
    end

    validIdx = idxInModel > 0;
    idxInModel = idxInModel(validIdx);
    kcatList.kcats = kcatList.kcats(validIdx);
    kcatList.kcatSource = kcatList.kcatSource(validIdx);
    kcatList.genes = kcatList.genes(validIdx);
    kcatList.substrates = kcatList.substrates(validIdx);
    kcatList.eccodes = kcatList.eccodes(validIdx);
    kcatList.wildcardLvl = kcatList.wildcardLvl(validIdx);
    kcatList.origin = kcatList.origin(validIdx);

    % ----------------- build a combined rank score ----------------------
    % Source priority rank (lower is better)
    srcRank = inf(numel(kcatList.kcats),1);
    for i = 1:numel(opts.priority)
        srcRank(strcmpi(kcatList.kcatSource, opts.priority{i})) = i;
    end

    % Fuzzy refinement: use origin (smaller is better) and wildcardLvl (smaller is better)
    isFuzzy = strcmpi(kcatList.kcatSource, 'brenda');
    origin  = getOr(kcatList, 'origin',      nan(numel(srcRank),1));
    wc      = getOr(kcatList, 'wildcardLvl', nan(numel(srcRank),1));

    origin2 = origin;   origin2(isnan(origin2)) = 999; % treat missing as worst
    wc2     = wc;       wc2(isnan(wc2))         = 999;

    % Linearized tuple ranking: source ≪ origin ≪ wildcard
    baseRank   = srcRank * 1e6;
    fuzzyRank  = zeros(size(baseRank));
    fuzzyRank(isFuzzy) = origin2(isFuzzy) * 1e3 + wc2(isFuzzy);
    comboRank  = baseRank + fuzzyRank;

    % Force CompleteMatch as best possible layer
    comboRank(strcmpi(kcatList.kcatSource,'CompleteMatch')) = -1;

    % --------------------- per-reaction selection -----------------------
    idxU = unique(idxInModel);
    selectedK   = zeros(numel(idxU),1);
    selectedSrc = cell(numel(idxU),1);

    for ii = 1:numel(idxU)
        ind  = idxU(ii);
        rows = find(idxInModel == ind);

        % Step 1: pick best layer by combined rank
        rmin  = min(comboRank(rows));
        rows1 = rows(comboRank(rows) == rmin);

        % Step 2: if some have genes, prefer those (traceability)
        if isfield(kcatList,'genes') && ~isempty(kcatList.genes)
            hasGene = ~cellfun(@isempty, kcatList.genes(rows1));
            if any(hasGene), rows1 = rows1(hasGene); end
        end

        % Step 3: tie-break inside the final set using the requested criterion
        vals = kcatList.kcats(rows1);
        switch lower(criteria)
            case 'max'
                target = max(vals);
                ties   = find(vals == target);
                if numel(ties) > 1
                    % deterministic tie-break: closest to median, then larger value, then lowest row index
                    medv = median(vals);
                    [~, order] = sortrows([abs(vals(ties) - medv), -vals(ties), rows1(ties)], [1 2 3]);
                    ties = ties(order(1));
                end
                pick = rows1(ties);
            case 'min'
                target = min(vals);
                pick   = rows1(find(vals == target, 1, 'first'));
            case 'median'
                m      = median(vals);
                [~, j] = min(abs(vals - m));
                target = m;
                pick   = rows1(j);
            case 'mean'
                m      = mean(vals);
                [~, j] = min(abs(vals - m));
                target = m;
                pick   = rows1(j);
            otherwise
                error('Invalid criteria: %s', criteria);
        end

        selectedK(ii)   = target;
        selectedSrc(ii) = kcatList.kcatSource(pick);
    end

    % -------------------------- write back ------------------------------
    switch lower(overwrite)
        case 'true'
            ec.kcat(idxU)   = selectedK;
            ec.source(idxU) = selectedSrc;
        case 'false'
            empty = (ec.kcat == 0);
            [idx2, which] = intersect(idxU, find(empty), 'stable');
            ec.kcat(idx2)   = selectedK(which);
            ec.source(idx2) = selectedSrc(which);
        case 'ifhigher'
            higher = ec.kcat(idxU) < selectedK;
            idx2   = idxU(higher);
            ec.kcat(idx2)   = selectedK(higher);
            ec.source(idx2) = selectedSrc(higher);
        otherwise
            error('Invalid overwrite flag: %s', overwrite);
    end

    model.enzymeConstraints = ec;
    rxnIdx = idxU;
end

% ----------------- helpers -----------------
function out = tern(c,a,b), if c, out=a; else, out=b; end, end
function v = getOr(S, f, defaultV)
    if isfield(S,f) && ~isempty(S.(f)), v = S.(f);
    else, v = defaultV; end
end
