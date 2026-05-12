function [rmse, expFlux, simFlux, matchResult] = rmsecal_13C(ecModel, C13Data, conditionIdx, bioRxn, cSource, C13ReactionMap, varargin)
% rmsecal_13C
%   Calculates RMSE between experimental and simulated 13C flux data.
%   Handles all three ecModel types: basic, isozyme, integrated.
%   Uses pre-built reaction mapping for efficiency.
%
% Input:
%   ecModel        - Enzyme-constrained model
%   C13Data        - 13C flux data (from load13CData)
%   conditionIdx   - Which condition to evaluate (index or name)
%   bioRxn         - Biomass reaction ID
%   cSource        - Carbon source exchange reaction ID
%   C13ReactionMap - Pre-built reaction mapping from buildC13ReactionMap
%   varargin       - Optional:
%       'verbose' - Print diagnostic messages (default: false)
%
% Output:
%   rmse        - RMSE value
%   expFlux     - Experimental flux values (matched reactions)
%   simFlux     - Simulated flux values (matched reactions)
%   matchResult - Reaction matching result

    p = inputParser;
    addOptional(p, 'verbose', false);
    parse(p, varargin{:});
    verbose = p.Results.verbose;

    % Parse condition
    if isnumeric(conditionIdx)
        condIdx = conditionIdx;
    elseif ischar(conditionIdx)
        condIdx = find(strcmp(C13Data.conditions, conditionIdx));
        if isempty(condIdx)
            error('Condition not found: %s', conditionIdx);
        end
    end

    % --- Validate that C13ReactionMap is provided ---
    if isempty(C13ReactionMap)
        error('C13ReactionMap is required. Build it using buildC13ReactionMap.');
    end

    % --- Step 1: Apply model constraints for this condition ---
    modelCond = applyConstraints(ecModel, C13Data, condIdx, bioRxn, cSource);

    % --- Step 2: Get experimental fluxes for this condition ---
    expFluxes = C13Data.fluxes(:, condIdx);
    validIdx = ~isnan(expFluxes);

    % Filter C13ReactionMap to valid indices
    nRxns = length(C13ReactionMap.reactions);
    validMask = validIdx(:)';  % Ensure row vector
    if length(validMask) < nRxns
        validMask(nRxns) = false;
    end

    if verbose
        fprintf('[rmsecal_13C] Condition %d: %d flux reactions\n', condIdx, sum(validMask));
    end

    % --- Step 3: Check if there are valid reactions to process ---
    nValidRxns = sum(validMask);
    if nValidRxns == 0
        warning('[rmsecal_13C] No valid 13C reactions for this condition');
        rmse = NaN;
        expFlux = [];
        simFlux = [];
        matchResult = struct();
        return;
    end

    % Get model type
    ecModelType = [];
    if isfield(ecModel, 'enzymeConstraints') && isfield(ecModel.enzymeConstraints, 'ecModeltype')
        ecModelType = ecModel.enzymeConstraints.ecModeltype;
    end
    isBasicModel = strcmp(ecModelType, 'basic');

    % Count matched reactions
    matchedMask = validMask;
    for i = 1:length(C13ReactionMap.modelIndices)
        if validMask(i) && isempty(C13ReactionMap.modelIndices{i})
            matchedMask(i) = false;
        end
    end
    nMatched = sum(matchedMask);

    if verbose
        fprintf('[rmsecal_13C] Condition %d: %d valid, %d matched\n', condIdx, nValidRxns, nMatched);
    end

    if nMatched == 0
        warning('[rmsecal_13C] No reactions matched between experiment and model');
        rmse = NaN;
        expFlux = [];
        simFlux = [];
        matchResult = struct();
        return;
    end

    % --- Step 4: Run FBA with pFBA ---
    % Find protein pool exchange reaction
    prot_idx = find(strcmp(modelCond.rxns, 'prot_pool_exchange'));

    modelCond.c(:) = 0;
    objIdx = find(strcmp(modelCond.rxns, bioRxn));
    if isempty(objIdx)
        objIdx = find(strcmp(modelCond.rxnNames, 'biomass'));
    end
    if isempty(objIdx)
        warning('[rmsecal_13C] Biomass reaction not found');
        rmse = 999;
        expFlux = [];
        simFlux = [];
        matchResult = struct();
        return;
    end

    % First FBA: maximize biomass
    modelCond.c(objIdx) = 1;
    sol = solveLP(modelCond);

    if isempty(sol.x) || sol.stat ~= 1
        % Silent fail - just return high RMSE for this condition
        rmse = 999;
        expFlux = [];
        simFlux = [];
        matchResult = struct();
        return;
    end
    
    % pFBA: fix biomass at 99% of optimum, minimize protein usage
    if ~isempty(prot_idx)
        modelCond.lb(objIdx) = sol.f * 0.99;
        modelCond.c(:) = 0;
        modelCond.c(prot_idx) = 1;  % Minimize protein pool
        sol = solveLP(modelCond);

        if isempty(sol.x) || sol.stat ~= 1
            if verbose
                warning('[rmsecal_13C] pFBA failed, using FBA solution');
            end
            % Keep original FBA solution
        end
    end

    % --- Step 5: Extract simulated fluxes and calculate carbon-weighted RMSE ---
    expCarbonFlux = zeros(nMatched, 1);
    simCarbonFlux = zeros(nMatched, 1);
    reactionNames = cell(nMatched, 1);
    resultIdx = 0;

    for i = 1:nRxns
        if ~validMask(i) || ~matchedMask(i)
            continue;
        end

        resultIdx = resultIdx + 1;

        % Get data from pre-built mapping
        reactionNames{resultIdx} = C13ReactionMap.reactions{i};
        modelIndices = C13ReactionMap.modelIndices{i};
        userDirections = C13ReactionMap.directions{i};  % User-provided directions from TSV

        % Get carbon count from pre-built mapping
        nCarbon = C13ReactionMap.nCarbon(i);

        % Experimental flux
        expNetFlux = expFluxes(i);

        % Calculate experimental carbon flux
        expCarbonFlux(resultIdx) = expNetFlux * nCarbon;

        % Aggregate simulated fluxes from matched reactions
        % User directions: {1} means use model forward direction, {-1} means reverse
        simNetFlux = 0;
        nReactions = length(modelIndices);

        for j = 1:nReactions
            idx = modelIndices(j);
            if idx <= length(sol.x)
                flux_j = sol.x(idx);
                % Get user-specified direction for this reaction
                if isempty(userDirections) || j > length(userDirections)
                    dir_j = 1;  % Default to forward
                else
                    dir_j = userDirections(j);
                end

                % Apply direction: 1 means use flux as-is, -1 means negate
                simNetFlux = simNetFlux + dir_j * flux_j;
            end
        end

        % Calculate simulated carbon flux
        simCarbonFlux(resultIdx) = simNetFlux * nCarbon;
    end

    expFlux = expCarbonFlux;
    simFlux = simCarbonFlux;

    % Calculate RMSE
    rmse = sqrt(mean((expFlux - simFlux).^2));

    if verbose
        fprintf('[rmsecal_13C] RMSE: %.4f (matched %d reactions)\n', rmse, nMatched);

        % Show detailed comparison for first few reactions
        nShow = min(5, length(reactionNames));
        fprintf('  Sample comparisons:\n');
        for j = 1:nShow
            rxn = reactionNames{j};
            if iscell(rxn)
                rxnStr = ['{', strjoin(rxn, ';'), '}'];
            else
                rxnStr = rxn;
            end
            fprintf('    %s: exp=%.3f, sim=%.3f (diff=%.3f)\n', ...
                rxnStr, expFlux(j), simFlux(j), abs(expFlux(j) - simFlux(j)));
        end
    end

    % Store match result for debugging
    matchResult = struct();
    matchResult.reactionNames = reactionNames;
    matchResult.expFlux = expFlux;
    matchResult.simFlux = simFlux;
    matchResult.nMatched = nMatched;
    matchResult.ecModelType = ecModelType;
end

%% ------------------- Helper Functions -------------------

function modelCond = applyConstraints(ecModel, C13Data, condIdx, bioRxn, cSource)
% applyConstraints
%   Applies environment constraints from 13C data to the model
%
%   Constraint rules:
%   - cSource (carbon source): always constrained
%   - Constraints with NaN values: set lb=-1000 for NaN, valid value for others
%   - Other constraints (biomass, etc.): NOT applied to model, only for RMSE
%
%   NaN values are not included in RMSE calculation

    modelCond = ecModel;
    modelCond.lb(strcmp(modelCond.rxns, cSource)) = 0;

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

            % Check if this constraint should be applied
            shouldConstrain = false;

            % 1. Always constrain cSource (carbon source)
            if strcmp(field, cSource)
                shouldConstrain = true;
            end

            % 2. Constrain if contains NaN (has partial valid data)
            if isfield(con, 'lb')
                lbVals = con.lb;
                if condIdx <= length(lbVals)
                    hasNaN = any(isnan(lbVals));
                    hasValid = any(~isnan(lbVals));
                    if hasNaN && hasValid
                        shouldConstrain = true;
                    end
                end
            end

            if ~shouldConstrain
                continue;  % Skip this constraint, only used for RMSE
            end

            % Apply the constraint
            if isfield(con, 'lb')
                lbVals = con.lb;
                if condIdx <= length(lbVals)
                    val = lbVals(condIdx);
                    if isnan(val)
                        % NaN: set to -1000 (essentially unconstrained)
                        modelCond.lb(rxnIdx) = -1000;
                    else
                        % Valid value: apply it
                        modelCond.lb(rxnIdx) = val;
                    end
                else
                    modelCond.lb(rxnIdx) = -1000;
                end
            end

            % Apply upper bound if available
            if isfield(con, 'ub')
                ubVals = con.ub;
                if condIdx <= length(ubVals)
                    val = ubVals(condIdx);
                    if isnan(val)
                        % NaN: set to 1000 (essentially unconstrained)
                        modelCond.ub(rxnIdx) = 1000;
                    else
                        % Valid value: apply it
                        modelCond.ub(rxnIdx) = val;
                    end
                else
                    modelCond.ub(rxnIdx) = 1000;
                end
            end
        end

        % Apply biomass objective if available (for RMSE calculation reference)
        if isfield(constraints, 'biomass')
            bioObj = constraints.biomass.objective;
            % Note: biomass constraint is NOT applied to lb here
            % Only used for RMSE calculation reference
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

