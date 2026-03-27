function equation = getReactionEquation(model, rxnIdx, varargin)
% getReactionEquation
%   Retrieves the reaction equation as a string for a given reaction.
%
% Input:
%   model       - Metabolic model structure
%   rxnIdx      - Index or name of the reaction
%   varargin    - Optional: 'metNames' to use metabolite names instead of IDs
%
% Output:
%   equation    - Reaction equation string (e.g., '2 GLC -> 2 G6P')
%
% Examples:
%   eq = getReactionEquation(model, 'PGI');
%   eq = getReactionEquation(model, 123);

    % Parse input
    p = inputParser;
    addOptional(p, 'useMetNames', false);
    parse(p, varargin{:});
    useMetNames = p.Results.useMetNames;

    % Get reaction index if name provided
    if ischar(rxnIdx)
        if isfield(model, 'rxns')
            rxnIdx = find(strcmp(model.rxns, rxnIdx));
        elseif isfield(model, 'rxnNames')
            rxnIdx = find(strcmp(model.rxnNames, rxnIdx));
        end
        if isempty(rxnIdx)
            error('Reaction not found: %s', rxnIdx);
        end
    end

    % Get metabolites involved in the reaction
    rxnS = model.S(:, rxnIdx);
    metIndices = find(rxnS ~= 0);

    if isempty(metIndices)
        equation = '';
        return;
    end

    % Build equation string
    reactants = {};
    products = {};

    for i = 1:length(metIndices)
        metIdx = metIndices(i);
        coeff = rxnS(metIdx);

        % Get metabolite identifier
        if useMetNames && isfield(model, 'metNames')
            metName = model.metNames{metIdx};
        else
            metName = model.mets{metIdx};
        end

        % Format coefficient (remove sign for display, will add later)
        if abs(coeff) == 1
            coeffStr = '';
        else
            coeffStr = num2str(abs(coeff), '%.2g');
        end

        metStr = [coeffStr metName];

        if coeff < 0
            reactants{end+1} = metStr;
        else
            products{end+1} = metStr;
        end
    end

    % Combine into equation
    if isempty(reactants)
        equation = ['-> ' strjoin(products, ' + ')];
    elseif isempty(products)
        equation = [strjoin(reactants, ' + ') + ' ->'];
    else
        equation = [strjoin(reactants, ' + ') + ' -> ' + strjoin(products, ' + ')];
    end
end
