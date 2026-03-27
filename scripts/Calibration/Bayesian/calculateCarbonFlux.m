function carbonFlux = calculateCarbonFlux(model, rxnIdx, flux, varargin)
% calculateCarbonFlux
%   Calculates carbon flux for a reaction by multiplying the flux value
%   by the total carbon atoms involved in the reaction.
%
% Input:
%   model       - Metabolic model structure
%   rxnIdx      - Index or name of the reaction
%   flux        - Flux value (e.g., from FBA solution)
%   varargin    - Optional:
%       'direction' - 'forward', 'backward', or 'net' (default: 'net')
%
% Output:
%   carbonFlux  - Carbon flux (flux * total carbon atoms)
%
% Note:
%   For reversible reactions, the carbon flux considers direction:
%   - Forward: carbon moving from reactants to products
%   - Backward: carbon moving from products to reactants
%   - Net: flux * sign (positive = forward, negative = backward)

    p = inputParser;
    addOptional(p, 'direction', 'net');
    parse(p, varargin{:});
    direction = p.Results.direction;

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

    % Get reaction stoichiometry
    rxnS = model.S(:, rxnIdx);

    % Calculate carbon atoms for each metabolite
    nCarbons = getMetaboliteCarbon(model);

    % Calculate total carbon in reactants and products
    reactantCarbons = sum(nCarbons(rxnS < 0) .* (-rxnS(rxnS < 0)));
    productCarbons = sum(nCarbons(rxnS > 0) .* rxnS(rxnS > 0));

    switch direction
        case 'forward'
            % Carbon flowing from reactants to products
            carbonFlux = flux * reactantCarbons;
        case 'backward'
            % Carbon flowing from products to reactants
            carbonFlux = flux * productCarbons;
        case 'net'
            % Net carbon flux (positive = forward, negative = backward)
            carbonFlux = flux * (reactantCarbons - productCarbons);
        otherwise
            error('Unknown direction: %s', direction);
    end
end


function nCarbons = getMetaboliteCarbon(model)
% getMetaboliteCarbon
%   Returns a vector of carbon atom counts for all metabolites in the model.

    nMets = length(model.mets);
    nCarbons = zeros(nMets, 1);

    % Check if metFormulas exists
    if ~isfield(model, 'metFormulas')
        warning('Model does not have metFormulas field');
        return;
    end

    % Parse formulas for carbon count
    for i = 1:nMets
        formula = model.metFormulas{i};
        if isempty(formula)
            nCarbons(i) = 0;
            continue;
        end

        % Extract carbon count (look for 'C' followed by number)
        carbonMatch = regexp(formula, 'C(\d*)', 'tokens');
        if ~isempty(carbonMatch)
            carbonStr = carbonMatch{1}{1};
            if isempty(carbonStr)
                nCarbons(i) = 1;
            else
                nCarbons(i) = str2double(carbonStr);
            end
        else
            nCarbons(i) = 0;
        end
    end
end


function carbonFluxes = calculateMultipleCarbonFluxes(model, fluxes, varargin)
% calculateMultipleCarbonFluxes
%   Calculates carbon fluxes for multiple reactions.
%
% Input:
%   model       - Metabolic model structure
%   fluxes      - Vector of flux values (or matrix for multiple solutions)
%   varargin    - Optional parameters passed to calculateCarbonFlux
%
% Output:
%   carbonFluxes - Carbon flux values (same size as fluxes input)

    % Get reaction indices (if flux is a vector matching rxns)
    if isvector(fluxes)
        nRxns = length(fluxes);
        carbonFluxes = zeros(size(fluxes));

        for i = 1:nRxns
            carbonFluxes(i) = calculateCarbonFlux(model, i, fluxes(i), varargin{:});
        end
    else
        error('Only vector flux input supported currently');
    end
end
