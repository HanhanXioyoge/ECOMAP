function [ecModel, finalKcats, initialKcats, rmseHistory, kcatHistory] = bayesianTuning(ecModel, useConstraint, useUnconstrained, useC13Flux, saveCurrentState, maxIterations, proc, numPerGeneration, rejectNum, initialKcatsInput, parameters)
    % 1. Initialization and Data Setup
    % -------------------------------------------------------------------------
    if nargin < 11 || isempty(parameters), parameters = ParameterManager.getParams(); end
    if nargin < 10 || isempty(initialKcatsInput), initialKcatsInput = ecModel.enzymeConstraints.kcat; end
    
    initialKcats = initialKcatsInput;
    nKcats  = length(initialKcats);
    logInit = log10(initialKcats);
    
    % Hard Bounds: +/- 3 orders of magnitude
    logLB_Hard = logInit - 3;
    logUB_Hard = logInit + 3;
    
    % Default parameters
    if nargin < 9 || isempty(rejectNum),       rejectNum = 0.2;        end
    if nargin < 8 || isempty(numPerGeneration),numPerGeneration = 160; end 
    if nargin < 7 || isempty(proc),             proc = parameters.PRESTO.ncpu; end
    if nargin < 6 || isempty(maxIterations),    maxIterations = 150;    end
    if nargin < 5 || isempty(saveCurrentState), saveCurrentState = true;end
    
    bioRxn = parameters.bioRxn; cSource = parameters.c_source;
    basePath = parameters.dataDir; analysisPath = parameters.outputDir;
    org_name = parameters.org_name;
    
    % [Fix] Pre-calculate Carbon Numbers
    if ~isfield(ecModel, 'excarbon')
        ecModel = addCarbonNum(ecModel, bioRxn);
    end
    
    % Determine data combination mode for state file naming
    % 4 possible combinations:
    %   1 = growth_con: useConstraint only
    %   2 = growth_con_uncon: useConstraint + useUnconstrained
    %   3 = growth_con_c13: useConstraint + useC13Flux
    %   4 = full_con_uncon_c13: all three data types
    if useConstraint && useUnconstrained && useC13Flux
        modeTag = 'full_con_uncon_c13';      % All three data types
    elseif useConstraint && useC13Flux
        modeTag = 'growth_con_c13';          % Constraint + 13C
    elseif useConstraint && useUnconstrained
        modeTag = 'growth_con_uncon';        % Constraint + Unconstrained
    else
        modeTag = 'growth_con';              % Constraint only (default)
    end
    stateFile = fullfile(analysisPath, ['bayesianState_LHS3_EDA_PCA_' modeTag '.mat']);
    
    % Data Loading
    growthData = [];  if useConstraint, growthData = readtable(fullfile(basePath,'BayesianGrowthRates.tsv'), 'FileType', 'text', 'ReadRowNames', true); end
    unconsData = [];  if useUnconstrained, unconsData = readtable(fullfile(basePath,'UnconstrainedMaxGrowth.tsv'), 'FileType', 'text', 'ReadRowNames', true); end

    % 13C Flux Data - Use dedicated loading function
    C13Fluxdata = [];
    if useC13Flux
        c13File = fullfile(basePath, '13CFluxdata.tsv');
        if exist(c13File, 'file')
            C13Fluxdata = load13CData(c13File, ecModel);
        else
            warning('13C flux data file not found: %s', c13File);
        end
    end
    
    % --- Robust Initialization ---
    kcats = initialKcats; 
    rmseHistory = zeros(0, 1); 
    kcatHistory = zeros(nKcats, 0); 
    sampledGeneration = 1;
    theta_100 = zeros(0, 1); 
    kcat_100  = zeros(nKcats, 0);
    
    % Stage Control: 1=LHS/EDA, 2=PCA
    currentStage = 1; 
    noImproveCount = 0;      
    lastBestRMSE = inf;
    lastWorstRMSE = inf;
    improvementThreshold = 0.001; 
    
    % 2. Resume Logic
    % -------------------------------------------------------------------------
    if saveCurrentState && exist(stateFile, 'file')
        fprintf('[Bayesian] State file found. Resuming...\n');
        vars = load(stateFile);
        if isfield(vars, 'currentStage'), currentStage = vars.currentStage; end
        load(stateFile, 'rmseHistory', 'kcatHistory', 'theta_100', 'kcat_100', 'sampledGeneration', 'initialKcats', 'lastBestRMSE', 'noImproveCount');
        
        % Ensure correct dimensions
        theta_100 = theta_100(:);
        sampledGeneration = sampledGeneration + 1;
        kcats = kcat_100(:, 1); 
        
        % Recalculate stats
        if length(theta_100) > 1
            BestRMSE = theta_100(1);
            WorstRMSE = theta_100(end);
        else
            BestRMSE = inf;
            WorstRMSE = inf;
        end
    else
        % Initial Run Check
        fprintf('[Init] Running initial baseline check...\n');
        [BestRMSE, ~, ~] = abc_max(ecModel, kcats, growthData, unconsData, C13Fluxdata, 1, 1, 1, bioRxn, cSource, [], org_name);
        WorstRMSE = inf;
        kcat_100 = kcats;
        theta_100 = BestRMSE;
    end
    
    % 3. Main Optimization Loop
    % -------------------------------------------------------------------------
    while sampledGeneration <= maxIterations && BestRMSE > rejectNum
        tic; 
        
        % Calculate Gap
        Gap_D = abs(WorstRMSE - BestRMSE);
        
        % --- Logic: Smart Switching based on Best AND Worst ---
        % Only switch if we are past the LHS phase (Gen > 3)
        if currentStage == 1 && sampledGeneration > 3
            % Criteria 1: Reached 90% of iterations -> Force PCA
            forceSwitch = sampledGeneration >= floor(maxIterations * 0.9);
            
            % Criteria 2: Stagnation Logic
            bestStagnated = (noImproveCount >= 5);
            worstStagnated = abs(lastWorstRMSE - WorstRMSE) < 1e-3;
            converged = (Gap_D < 1); 
            
            if forceSwitch
                currentStage = 2;
                noImproveCount = 0;
                fprintf('  [Switch] 90%% iterations done. Forcing Stage 2 (PCA)...\n');
            elseif bestStagnated && sampledGeneration > 5
                if converged || worstStagnated
                    currentStage = 2;
                    noImproveCount = 0;
                    fprintf('  [Switch] Stagnation detected (Best & Worst stuck). Switching to Stage 2 (PCA)...\n');
                else
                    fprintf('  [Info] Best stuck, but Worst is improving. Continuing EDA...\n');
                    noImproveCount = 0; 
                end
            end
        end
        
        % Method Name Logging
        if sampledGeneration <= 3
            methodName = 'LHS-Init';
        elseif currentStage == 1
            methodName = 'Normal-EDA'; 
        else
            methodName = 'PCA-Corr'; 
        end
        
        fprintf('\n[Bayesian-%s] Gen %d/%d [%s]\n', modeTag, sampledGeneration, maxIterations, methodName);
        fprintf('  Best: %.4f | Worst: %.4f | Gap: %.4f\n', BestRMSE, WorstRMSE, Gap_D);
        
        % --- SAMPLING PHASE ---
        if sampledGeneration <= 3
            %% PHASE 1: Triple-Layer LHS Initialization
            % Gen 1: +/- 3 range (1024 samples) - Global Search
            % Gen 2: +/- 2 range (1024 samples) - Regional Search
            % Gen 3: +/- 1 range (1024 samples) - Local Search
            
            if sampledGeneration == 1
                rangeVal = 2.5;
                fprintf('  [Gen 1] LHS: Global Search (+/- 2.5 range, 1024 samples)...\n');
            elseif sampledGeneration == 2
                rangeVal = 2;
                fprintf('  [Gen 2] LHS: Regional Search (+/- 2 range, 1024 samples)...\n');
            else
                rangeVal = 1.5;
                fprintf('  [Gen 3] LHS: Local Search (+/- 1.5 range, 1024 samples)...\n');
            end
            
            currentSize = 1024;
            logLB_curr = logInit - rangeVal; 
            logUB_curr = logInit + rangeVal;
            
            lhsRaw = lhsdesign(currentSize, nKcats);
            % Linear scaling to current bounds
            kcatRandomAll = 10.^(logLB_curr + (logUB_curr - logLB_curr) .* lhsRaw');
            clear lhsRaw;
            
        elseif currentStage == 1
            %% PHASE 2: EDA (Normal Distribution based on Elite Population)
            fprintf('  [Stage 1] EDA: Normal Sampling based on Elite Population (N=%d)...\n', numPerGeneration);
            
            kcatRandomAll = zeros(nKcats, numPerGeneration);
            
            logElites = log10(kcat_100);
            mu_log = mean(logElites, 2);     
            sigma_log = std(logElites, 0, 2); 
            
            % Minimal variance protection
            sigma_log = max(sigma_log, 0.05);
            
            % Generate independent samples
            kcatRandomAll = 10.^(mu_log + randn(nKcats, numPerGeneration) .* sigma_log);
            kcatRandomAll = max(min(kcatRandomAll, 10.^logUB_Hard), 10.^logLB_Hard);
            
        else
            %% PHASE 3: PCA (Correlated Search)
            fprintf('  [Stage 2] PCA: Searching for metabolic correlations...\n');
            
            logElite = log10(kcat_100);
            bestLogKcat = logElite(:, 1);
            
            [coeff, ~, latent] = pca(logElite', 'Economy', true);
            numPC = find(cumsum(latent)/sum(latent) > 0.95, 1);
            numPC = max(30, min(numPC, 50)); 
            
            lStd = sqrt(latent(1:numPC));
            lStd = max(lStd, 0.05); 
            
            shrinkFactor = 0.5; 
            coeffRed = coeff(:, 1:numPC);
            
            scores = randn(numPerGeneration, numPC) .* (lStd' * shrinkFactor);
            logNew = repmat(bestLogKcat, 1, numPerGeneration) + coeffRed * scores';
            
            logNew = max(min(logNew, logUB_Hard), logLB_Hard);
            kcatRandomAll = 10.^logNew;
        end
        
        % 4. Evaluation and Elite Selection (Robust)
        % -------------------------------------------------------------------------
        batchCount = ceil(size(kcatRandomAll, 2) / proc);
        kcatBatches = cell(proc, 1);
        for j = 1:proc
            sIdx = (j-1)*batchCount + 1; eIdx = min(j*batchCount, size(kcatRandomAll, 2));
            if sIdx <= eIdx, kcatBatches{j} = kcatRandomAll(:, sIdx:eIdx); end
        end
        
        newRmseCells = cell(proc, 1);
        parfor j = 1:proc
            if ~isempty(kcatBatches{j})
                [newRmseCells{j}, ~, ~] = abc_max(ecModel, kcatBatches{j}, growthData, unconsData, C13Fluxdata, 1, size(kcatBatches{j}, 2), 1, bioRxn, cSource, [], org_name);
            end
        end
        
        % [Fix] Correct Merging Logic (Horizontal Concatenation)
        newRmse = [newRmseCells{:}]'; 
        
        allTheta = [newRmse; theta_100];
        allKcats = [kcatRandomAll, kcat_100];
        
        [~, sortedIdx] = sort(allTheta, 'ascend');
        nKeep = min(100, length(allTheta));
        theta_100 = allTheta(sortedIdx(1:nKeep));
        kcat_100  = allKcats(:, sortedIdx(1:nKeep));
        
        % Update Stats
        kcats = kcat_100(:, 1);
        BestRMSE = theta_100(1);        
        WorstRMSE = theta_100(end);     
        
        % Check for Improvement (Based on Best RMSE)
        if abs(lastBestRMSE - BestRMSE) < improvementThreshold
            noImproveCount = noImproveCount + 1;
        else
            noImproveCount = 0; 
        end
        
        % History Tracking
        lastBestRMSE = BestRMSE;
        lastWorstRMSE = WorstRMSE; 
        
        rmseHistory = [rmseHistory; BestRMSE];
        kcatHistory = [kcatHistory, kcats(:)];
        
        if saveCurrentState
            save(stateFile, 'rmseHistory', 'kcatHistory', 'theta_100', 'kcat_100', 'sampledGeneration', 'initialKcats', 'lastBestRMSE', 'noImproveCount', 'currentStage');
        end
        
        genTime = toc;
        fprintf('  [Result] Time: %.1fs | Best: %.4f | Change: %.2e\n', ...
                genTime, BestRMSE, abs(lastBestRMSE - BestRMSE));
        
        % --- Final Termination Logic ---
        if currentStage == 2 && noImproveCount >= 5
            fprintf('\n  [Termination] Final convergence in PCA stage. No improvement for 5 generations.\n');
            break; 
        end
        
        sampledGeneration = sampledGeneration + 1;
        clear newRmseCells allTheta allKcats kcatRandomAll kcatBatches;
    end
    
    finalKcats = kcats;
    ecModel.enzymeConstraints.kcat = finalKcats;
end
function y = ternary(cond, a, b), if cond, y = a; else, y = b; end, end
