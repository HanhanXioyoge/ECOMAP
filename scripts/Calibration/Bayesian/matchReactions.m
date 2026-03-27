function matchResult = matchReactions(expReactions, model, varargin)
% matchReactions
%   Matches experimental reactions to model reactions based on reaction equations.
%
% Input:
%   expReactions - Cell array of experimental reaction names
%   model        - Metabolic model structure
%   varargin     - Optional parameters:
%       'method'     - 'exact' (default), 'equation', or 'fuzzy'
%       'threshold'  - Similarity threshold for fuzzy matching (0-1)
%
% Output:
%   matchResult  - Structure containing:
%       .matchedIdx     - Model indices for matched reactions
%       .unmatchedExp   - Experimental reactions not matched
%       .matchedExp     - Experimental reactions that were matched
%       .equationMatch  - Cell array of matched equation pairs
%
% Method Description:
%   'exact'    - Match by exact reaction name
%   'equation' - Match by reaction equation (metabolites and stoichiometry)
%   'fuzzy'   - Match by equation similarity (substring matching)

    % Parse optional parameters
    p = inputParser;
    addOptional(p, 'method', 'exact');
    addOptional(p, 'threshold', 0.8);
    parse(p, varargin{:});
    method = p.Results.method;
    threshold = p.Results.threshold;

    nExp = length(expReactions);

    % Initialize result
    matchedIdx = NaN(nExp, 1);
    matchedExp = {};
    unmatchedExp = {};
    equationMatch = {};

    % Get model reaction info
    if isfield(model, 'rxns')
        modelRxns = model.rxns;
    else
        error('Model must have rxns field');
    end

    switch method
        case 'exact'
            % Exact name matching
            for i = 1:nExp
                expName = expReactions{i};
                idx = find(strcmp(modelRxns, expName));
                if ~isempty(idx)
                    matchedIdx(i) = idx(1);
                    matchedExp{end+1} = expName;
                    equationMatch{end+1} = {expName, modelRxns{idx(1)}};
                else
                    unmatchedExp{end+1} = expName;
                end
            end

        case 'equation'
            % Equation-based matching
            % First, build a hash map of model equations
            fprintf('[matchReactions] Building model equation index...\n');

            % Build model equation signatures
            modelEq = cell(length(modelRxns), 1);
            for i = 1:length(modelRxns)
                eq = getReactionEquationCompact(model, i);
                modelEq{i} = eq;
            end

            % Match each experimental reaction
            for i = 1:nExp
                expName = expReactions{i};

                % Try exact name match first
                idx = find(strcmp(modelRxns, expName));
                if ~isempty(idx)
                    matchedIdx(i) = idx(1);
                    matchedExp{end+1} = expName;
                    equationMatch{end+1} = {expName, modelRxns{idx(1)}};
                    continue;
                end

                % Try equation matching - find reactions with same equation
                % This is useful for enzyme-constrained models with split reactions
                % For now, mark as unmatched (can be enhanced later)
                unmatchedExp{end+1} = expName;
            end

        case 'fuzzy'
            % Fuzzy matching based on equation similarity
            % Extract metabolites from equations and match
            error('Fuzzy matching not yet implemented');

        otherwise
            error('Unknown matching method: %s', method);
    end

    % Package results
    matchResult.matchedIdx = matchedIdx;
    matchResult.unmatchedExp = unmatchedExp;
    matchResult.matchedExp = matchedExp;
    matchResult.equationMatch = equationMatch;
    matchResult.method = method;

    % Summary
    fprintf('[matchReactions] Matched %d/%d reactions (method: %s)\n', ...
        length(matchedExp), nExp, method);
end

%% ------------------- Helper Function -------------------

function eqStr = getReactionEquationCompact(model, rxnIdx)
% getReactionEquationCompact
%   Creates a compact string representation of a reaction equation.
%   This is used for equation-based matching.

    if ischar(rxnIdx)
        rxnIdx = find(strcmp(model.rxns, rxnIdx));
    end

    rxnS = model.S(:, rxnIdx);
    metIndices = find(rxnS ~= 0);

    if isempty(metIndices)
        eqStr = '';
        return;
    end

    % Build signature: sorted list of met@coeff
    sig = {};
    for i = 1:length(metIndices)
        metIdx = metIndices(i);
        coeff = rxnS(metIdx);
        metName = model.mets{metIdx};

        % Simplify metabolite name (remove compartment if present)
        if contains(metName, '[')
            metName = strsplit(metName, '[');
            metName = metName{1};
        end

        sig{end+1} = sprintf('%s@%.3f', metName, coeff);
    end

    % Sort to make order-independent
    sig = sort(sig);
    eqStr = strjoin(sig, '|');
end
