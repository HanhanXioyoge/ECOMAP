function [rmse, expFlux, simFlux, matchResult] = rmsecal_13C(ecModel, C13Data, conditionIdx, bioRxn, cSource, varargin)
% rmsecal_13C
%   Calculates RMSE between experimental and simulated 13C flux data.
%
% Input:
%   ecModel      - Enzyme-constrained model
%   C13Data     - 13C flux data (from load13CData)
%   conditionIdx - Which condition to evaluate (index or name)
%   bioRxn      - Biomass reaction ID
%   cSource     - Carbon source exchange reaction ID
%   varargin    - Optional:
%       'useCarbonFlux' - Use carbon-weighted flux (default: true)
%       'normalize'     - Normalize by substrate uptake (default: true)
%
% Output:
%   rmse        - RMSE value
%   expFlux     - Experimental flux values (matched reactions)
%   simFlux     - Simulated flux values (matched reactions)
%   matchResult - Reaction matching result

    p = inputParser;
    addOptional(p, 'useCarbonFlux', true);
    addOptional(p, 'normalize', true);
    parse(p, varargin{:});
    useCarbonFlux = p.Results.useCarbonFlux;
    normalize = p.Results.normalize;

    % Parse condition
    if isnumeric(conditionIdx)
        condIdx = conditionIdx;
    elseif ischar(conditionIdx)
        condIdx = find(strcmp(C13Data.conditions, conditionIdx));
        if isempty(condIdx)
            error('Condition not found: %s', conditionIdx);
        end
    end

    % --- Step 1: Apply model constraints for this condition ---
    modelCond = applyConstraints(ecModel, C13Data, condIdx, bioRxn, cSource);

    % --- Step 2: Get experimental fluxes for this condition ---
    expFluxes = C13Data.fluxes(:, condIdx);
    validIdx = ~isnan(expFluxes);
    expReactions = C13Data.reactions(validIdx);
    expFlux = expFluxes(validIdx);

    fprintf('[rmsecal_13C] Condition %d: %d reactions\n', condIdx, length(expReactions));

    % --- Step 3: Match reactions to model ---
    % Use pre-computed matching if available
    if isfield(C13Data, 'matchedIdx')
        matchedModelIdx = C13Data.matchedIdx(validIdx);
        matchedModelIdx = matchedModelIdx(matchedModelIdx > 0);
    else
        % Fallback: match by name
        matchResult_temp = matchReactions(expReactions, modelCond, 'method', 'exact');
        matchedModelIdx = matchResult_temp.matchedIdx(~isnan(matchResult_temp.matchedIdx));
    end

    if isempty(matchedModelIdx)
        warning('No reactions matched between experiment and model');
        rmse = NaN;
        simFlux = [];
        matchResult = struct();
        return;
    end

    % --- Step 4: Run FBA ---
    modelCond.c(:) = 0;
    objIdx = find(strcmp(modelCond.rxns, bioRxn));
    if isempty(objIdx)
        objIdx = find(strcmp(modelCond.rxnNames, 'biomass'));
    end
    if isempty(objIdx)
        warning('Biomass reaction not found');
        rmse = 999;
        simFlux = [];
        matchResult = struct();
        return;
    end
    modelCond.c(objIdx) = 1;

    sol = solveLP(modelCond);

    if isempty(sol.x) || sol.stat ~= 1
        warning('FBA solution failed');
        rmse = 999;
        simFlux = NaN(length(matchedModelIdx), 1);
        return;
    end

    % --- Step 5: Extract simulated fluxes and calculate RMSE ---
    % Get the experimental data for matched reactions
    matchedExpFlux = expFlux(validIdx);
    matchedMets = C13Data.metabolites(validIdx);
    matchedCoefs = C13Data.coefficients(validIdx);

    % Calculate net carbon flux RMSE
    expCarbonFlux = zeros(length(matchedExpFlux), 1);
    simCarbonFlux = zeros(length(matchedExpFlux), 1);

    for i = 1:length(matchedExpFlux)
        % Get reaction stoichiometry
        mets = matchedMets{i};
        coefs = matchedCoefs{i};

        % Calculate total carbon for reactants and products
        reactantC = 0;
        productC = 0;

        for j = 1:length(mets)
            metId = mets{j};
            coef = coefs(j);

            % Find carbon number for this metabolite
            metIdx = find(strcmp(modelCond.mets, metId));
            if ~isempty(metIdx)
                formula = modelCond.metFormulas{metIdx};
                cNum = extractCarbonFromFormula(formula);
            else
                cNum = 1; % Default
            end

            if coef > 0
                reactantC = reactantC + coef * cNum;
            else
                productC = productC + abs(coef) * cNum;
            end
        end

        % Calculate carbon-weighted flux
        if useCarbonFlux
            % Experimental flux (normalized, multiply by carbon)
            expCarbonFlux(i) = matchedExpFlux(i) * reactantC;

            % Simulated flux - need to get reaction flux from FBA solution
            % For now, use the reaction flux directly
            if matchedModelIdx(i) <= length(sol.x)
                rxnFlux = sol.x(matchedModelIdx(i));
                if rxnFlux >= 0
                    simCarbonFlux(i) = rxnFlux * reactantC;
                else
                    simCarbonFlux(i) = rxnFlux * productC;
                end
            else
                simCarbonFlux(i) = 0;
            end
        else
            expCarbonFlux(i) = matchedExpFlux(i);
            if matchedModelIdx(i) <= length(sol.x)
                simCarbonFlux(i) = sol.x(matchedModelIdx(i));
            else
                simCarbonFlux(i) = 0;
            end
        end
    end

    expFlux = expCarbonFlux;
    simFlux = simCarbonFlux;

    % Calculate RMSE
    rmse = sqrt(mean((expFlux - simFlux).^2));

    fprintf('[rmsecal_13C] RMSE: %.4f (matched %d reactions)\n', rmse, length(expFlux));

    matchResult = struct();
end

%% ------------------- Helper Functions -------------------

function modelCond = applyConstraints(ecModel, C13Data, condIdx, bioRxn, cSource)
% applyConstraints
%   Applies environment constraints from 13C data to the model

    modelCond = ecModel;

    % Apply constraints if available
    if isfield(C13Data, 'constraints') && ~isempty(fieldnames(C13Data.constraints))
        constraints = C13Data.constraints;
        constraintFields = fieldnames(constraints);

        for i = 1:length(constraintFields)
            field = constraintFields{i};
            con = constraints.(field);

            % Find reaction index
            rxnIdx = find(strcmp(modelCond.rxns, field));

            if isempty(rxnIdx)
                continue;
            end

            % Apply bounds
            if isfield(con, 'lb')
                lbVals = con.lb;
                if condIdx <= length(lbVals) && ~isnan(lbVals(condIdx))
                    modelCond.lb(rxnIdx) = lbVals(condIdx);
                end
            end

            if isfield(con, 'ub')
                ubVals = con.ub;
                if condIdx <= length(ubVals) && ~isnan(ubVals(condIdx))
                    modelCond.ub(rxnIdx) = ubVals(condIdx);
                end
            end
        end

        % Apply biomass objective if available
        if isfield(constraints, 'biomass')
            bioObj = constraints.biomass.objective;
            if condIdx <= length(bioObj) && ~isnan(bioObj(condIdx))
                bioRxnIdx = find(strcmp(modelCond.rxns, bioRxn));
                if ~isempty(bioRxnIdx)
                    modelCond.lb(bioRxnIdx) = bioObj(condIdx);
                end
            end
        end
    else
        % Default: set up aerobic glucose condition
        % Block other carbon sources
        modelCond.lb(strcmp(modelCond.rxns, cSource)) = 0;
        % Set glucose uptake
        glcIdx = find(strcmp(modelCond.rxns, 'EX_glc__D_e'));
        if ~isempty(glcIdx)
            modelCond.lb(glcIdx) = -10; % Default
        end
    end
end


function nC = extractCarbonFromFormula(formula)
% extractCarbonFromFormula
%   Extracts the number of carbon atoms from a molecular formula.

    if isempty(formula)
        nC = 0;
        return;
    end

    carbonMatch = regexp(formula, 'C(\d*)', 'tokens');
    if ~isempty(carbonMatch)
        carbonStr = carbonMatch{1}{1};
        if isempty(carbonStr)
            nC = 1;
        else
            nC = str2double(carbonStr);
        end
    else
        nC = 0;
    end
end
