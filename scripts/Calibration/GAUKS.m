function [ecModel_aerobic, ecModel_anaerobic, summaryTbl] = GAUKS(ecModel, AnaerobicModel, parameters)
% GAUKS: Global Aerobic/Anaerobic Uptake Kinetic Calibration
% This function calibrates enzyme-constrained models by setting Umin values
% for pre-structured exchange reactions (from applySluiceStructure).
%
% IMPORTANT: This function assumes the sluice structure has already been
% applied using applySluiceStructure(). It only probes Umin and calibrates.
%
% Usage:
%   1. Apply structure BEFORE calibration:
%      [model_sluice, config] = applySluiceStructure(model, ex_rxn_list);
%
%   2. Run Bayesian calibration:
%      finalKcats = bayesianTuning(model_sluice, ...);
%
%   3. Run GAUKS to calibrate Umin:
%      [model_aerobic, model_anaerobic, summary] = GAUKS(model_sluice, true, params);
%
% Note: Structure is fixed during calibration, so multiple kcat sets can be
%       applied to the same structured model.

%% 1. Initialization and Parameter Setup
    if nargin < 3 || isempty(parameters), parameters = ParameterManager.getParams(); end
    if nargin < 2 || isempty(AnaerobicModel), AnaerobicModel = false; end

    bioRxn      = parameters.bioRxn;
    basePath    = parameters.dataDir;
    prot_pool   = 'prot_pool';
    c_source    = parameters.c_source;
    org_name    = parameters.org_name;

    % Load experimental data
    T = readtable(fullfile(basePath, 'UnconstrainedMaxGrowth.tsv'), 'FileType','text','ReadRowNames',true);

    % Initialize Summary Table with new result columns
    summaryTbl = T;
    summaryTbl.Umin           = zeros(height(T), 1);
    summaryTbl.ProteomicCost  = zeros(height(T), 1); % This is the 'xi' value (proteomic cost)
    summaryTbl.Error_Percent  = zeros(height(T), 1);

    ex_o2              = T.Properties.VariableNames{9};
    exchange_reactions = T{:,1};
    TargetMu           = T{:,3};
    Oxavail            = T{:,11};

    % Clear carbon source
    idxC = strcmp(ecModel.rxns, c_source);
    if any(idxC), ecModel.lb(idxC) = 0; end

    uniqueExRxns = unique(exchange_reactions);
    ecModel_aerobic = ecModel;
    ecModel_anaerobic = [];

    if AnaerobicModel
        initial_anaerobic = anaerobicModel(ecModel, org_name, ex_o2);
        ecModel_anaerobic = initial_anaerobic;
    end

    %% Pre-check: Verify all conditions have valid solutions
    fprintf('[GAUKS] Pre-checking all conditions for feasibility...\n');
    failedConditions = {};

    for i = 1:numel(uniqueExRxns)
        ex_rxn = uniqueExRxns{i};

        % Check aerobic
        idxAeroRows = find(strcmp(Oxavail, 'aerobic') & strcmp(exchange_reactions, ex_rxn));
        if ~isempty(idxAeroRows)
            currentMu = max(TargetMu(idxAeroRows));
            if ~checkFeasibility(ecModel, ex_rxn, currentMu, 'aerobic', bioRxn, ex_o2)
                failedConditions{end+1} = sprintf('%s (aerobic, Mu=%.4f)', ex_rxn, currentMu); %#ok<AGROW>
            end
        end

        % Check anaerobic
        if AnaerobicModel
            idxAnaRows = find(strcmp(Oxavail, 'anaerobic') & strcmp(exchange_reactions, ex_rxn));
            if ~isempty(idxAnaRows)
                mu_anaero = max(TargetMu(idxAnaRows));
                if ~checkFeasibility(initial_anaerobic, ex_rxn, mu_anaero, 'anaerobic', bioRxn, ex_o2)
                    failedConditions{end+1} = sprintf('%s (anaerobic, Mu=%.4f)', ex_rxn, mu_anaero); %#ok<AGROW>
                end
            end
        end
    end

    % If any failed, throw error
    if ~isempty(failedConditions)
        error('GAUKS:InfeasibleConditions', ...
            'The following conditions are infeasible:\n%s\n', ...
            strjoin(failedConditions, '\n'));
    end
    fprintf('[GAUKS] All conditions are feasible. Proceeding...\n');

    %% 2. Phase 1: Probing Umin (Structure already applied!)
    % Note: Structure is assumed to be pre-applied via applySluiceStructure()
    for i = 1:numel(uniqueExRxns)
        ex_rxn = uniqueExRxns{i};

        % --- Aerobic Branch ---
        idxAeroRows = find(strcmp(Oxavail, 'aerobic') & strcmp(exchange_reactions, ex_rxn));
        if ~isempty(idxAeroRows)
            fprintf('[GAUKS] Probing Umin: %s (Aerobic)\n', ex_rxn);
            currentMu = max(TargetMu(idxAeroRows));
            umin_aero = probeUmin(ecModel, ex_rxn, currentMu, 'aerobic', bioRxn, ex_o2);
            % Record Umin
            summaryTbl.Umin(idxAeroRows) = umin_aero;
            % Apply Umin to the pre-structured model
            ecModel_aerobic = setUmin(ecModel_aerobic, ex_rxn, umin_aero, prot_pool);
        end

        % --- Anaerobic Branch ---
        if AnaerobicModel
            idxAnaRows = find(strcmp(Oxavail, 'anaerobic') & strcmp(exchange_reactions, ex_rxn));
            if ~isempty(idxAnaRows)
                fprintf('[GAUKS] Probing Umin: %s (Anaerobic)\n', ex_rxn);
                mu_anaero = max(TargetMu(idxAnaRows));
                umin_ana = probeUmin(initial_anaerobic, ex_rxn, mu_anaero, 'anaerobic', bioRxn, ex_o2);
                summaryTbl.Umin(idxAnaRows) = umin_ana;
                % Apply Umin to the pre-structured model
                ecModel_anaerobic = setUmin(ecModel_anaerobic, ex_rxn, umin_ana, prot_pool);
            end
        end
    end

    %% 3. Phase 3: Bisection Calibration
    fprintf('[GAUKS] Starting Bisection Calibration...\n');
    
    % Calibrate Aerobic Model and update Summary Table
    [ecModel_aerobic, summaryTbl] = calibrateModelMode(ecModel_aerobic, summaryTbl, 'aerobic', bioRxn);
    
    % Calibrate Anaerobic Model and update Summary Table
    if AnaerobicModel
        [ecModel_anaerobic, summaryTbl] = calibrateModelMode(ecModel_anaerobic, summaryTbl, 'anaerobic', bioRxn);
    end

    fprintf('[GAUKS] Calibration Complete.\n');
end

%% --- Helper 1: Double Optimization for Umin ---
function feasible = checkFeasibility(model, ex_rxn, targetMu, mode, bioRxn, ex_o2)
% CHECKFEASIBILITY
%   Check if a condition is feasible (has solution) before running GAUKS.
%
% Input:
%   model     - Enzyme-constrained model
%   ex_rxn   - Exchange reaction name
%   targetMu  - Target growth rate
%   mode      - 'aerobic' or 'anaerobic'
%   bioRxn   - Biomass reaction name
%   ex_o2    - Oxygen exchange reaction name
%
% Output:
%   feasible  - true if feasible, false otherwise

    % Apply oxygen constraint
    testModel = model;
    if strcmp(mode, 'anaerobic')
        testModel.lb(strcmp(testModel.rxns, ex_o2)) = 0;
    end

    % Fix growth rate
    idxBio = strcmp(testModel.rxns, bioRxn);
    testModel.lb(idxBio) = targetMu;
    testModel.lb(strcmp(testModel.rxns, ex_rxn)) = -1000;

    sol = solveLP(testModel);

    feasible = ~isempty(sol.x) && sol.stat == 1;
end

function umin = probeUmin(model, ex_rxn, targetMu, mode, bioRxn, ex_o2)
    % Enforce O2 constraints for anaerobic probing
    if strcmp(mode, 'anaerobic'), model.lb(strcmp(model.rxns, ex_o2)) = 0; end

    % Fix growth rate with high precision
    idxBio = strcmp(model.rxns, bioRxn);
    model.lb(idxBio) = targetMu * 0.999;
    model.ub(idxBio) = targetMu * 1.001;

    % Step 1: Minimize substrate absorption (Maximize exchange flux)
    model.c(:) = 0;
    model.lb(strcmp(model.rxns, ex_rxn)) = -1000;

    model.c(strcmp(model.rxns, ex_rxn)) = 1;
    sol1 = solveLP(model);
    
    % Relax protein pool if metabolic constraints prevent reaching TargetMu
    if isempty(sol1.x) || sol1.stat ~= 1
        pIdx = strcmp(model.rxns, 'prot_pool_exchange');
        if any(pIdx), model.lb(pIdx) = -1000; model.ub(pIdx) = 1000; end
        sol1 = solveLP(model);
    end

    if isempty(sol1.x) || sol1.stat ~= 1, umin = 1e-5; return; end
    umin_raw = abs(sol1.x(strcmp(model.rxns, ex_rxn)));

    % Step 2: Minimize protein usage at Step 1's uptake rate
    idxEx = strcmp(model.rxns, ex_rxn);
    model.lb(idxEx) = -umin_raw * 1.001;
    model.ub(idxEx) = -umin_raw * 0.999;
    model.c(:) = 0;
    pIdx = strcmp(model.rxns, 'prot_pool_exchange');
    if any(pIdx), model.c(pIdx) = 1; else, error('Missing protein pool exchange.'); end

    sol2 = solveLP(model);
    if ~isempty(sol2.x) && sol2.stat == 1, umin = abs(sol2.x(idxEx)); else, umin = umin_raw; end
end

%% --- Helper 2: Calibration Loop ---
function [model, T] = calibrateModelMode(model, T, mode, bioRxn)
    idxMode = find(strcmp(T{:,11}, mode)); % Using column 'Condition' (assume 11th col)
    relevantEx = unique(T{:,1}(idxMode));
    
    for i = 1:numel(relevantEx)
        ex_rxn = relevantEx{i};
        idxRows = idxMode(strcmp(T{:,1}(idxMode), ex_rxn));
        targetMu = max(T{:,3}(idxRows));
        
        fprintf('[GAUKS] Calibrating %s for %s mode (TargetMu: %.3f)...\n', ex_rxn, mode, targetMu);
        
        [model, bestXi, finalMu] = bisectionXi(model, ex_rxn, targetMu, bioRxn);
        
        % Populate Summary Table
        T.ProteomicCost(idxRows) = bestXi;
        T.Error_Percent(idxRows) = abs(finalMu - targetMu) / targetMu * 100;
    end
end

%% --- Helper 4: Bisection Core ---
function [model, xi, predictedMu] = bisectionXi(model, ex_rxn, targetMu, bioRxn)
    % BISECTIONXI: Bisection method to calibrate xi (proteomic cost)
    %   xi is the coefficient of prot_pool in the extended branch
    tol = 1e-4; maxIter = 100; low = 0; high = 10;
    extIdx  = find(strcmp(model.rxns, [ex_rxn '_extended']));
    protRow = find(strcmp(model.mets, 'prot_pool'));
    bioIdx  = strcmp(model.rxns, bioRxn);

    % Simulate excess substrate
    originalLB = model.lb(strcmp(model.rxns, ex_rxn));
    model.lb(strcmp(model.rxns, ex_rxn)) = -1000;

    xi = low; % Initialize
    predictedMu = 0;

    for iter = 1:maxIter
        xi = (low + high) / 2;
        model.S(protRow, extIdx) = xi;
        sol = solveLP(model);

        if ~isempty(sol.x) && sol.stat == 1
            predictedMu = sol.x(bioIdx);
        else
            predictedMu = 0;
        end

        if abs(predictedMu - targetMu) < tol
            fprintf('       Converged: Iter %d, Xi = %.6f, Err = %.2f%%\n', ...
                    iter, xi, abs(predictedMu-targetMu)/targetMu*100);
            break;
        end

        if predictedMu < targetMu
            low = xi; 
        else
            high = xi; 
        end
    end

    % Restore Switch
    model.lb(strcmp(model.rxns, ex_rxn)) = originalLB;
end