function result = algOptforce(model, biomassRxn, targetRxn, CsourceRxn, outputFile,parameters)
% algOptforce  Find forced flux changes using optForce (Mendoza et al.).
%   result = algOptforce(model, biomassRxn, targetRxn, CsourceRxn, outputFile)
%
%   Mirrors the optForce procedure:
%     1. FBA on WT (max biomass) and target (max target)
%     2. FVA on WT and mutant (constrained) strains
%     3. findMustU / findMustL (first-order must sets)
%     4. findMustUU / findMustLL / findMustUL (second-order must sets)
%     5. optForce (forced intervention sets, K=1 then K=2)
%
%   All mustU / mustL / mustUU / mustLL / mustUL sets are stored in result.targets.
%   Re-run overwrites previous optforce.* files in the workspace subfolder.
%   Per-set interventions are grouped together in the CSV; sets are separated
%   by a blank line.
%
%   Inputs:
%     model       struct, COBRA-compatible model
%     biomassRxn  char, biomass reaction ID
%     targetRxn   char, target reaction ID
%     CsourceRxn  char, carbon source uptake reaction ID
%     outputFile  char, output CSV path (empty -> stdout; default -> workspace)
%
%   Output:
%     result - struct (algFseof-style 7 fields):
%         config     - actual params used
%         biomassRxn - biomass reaction ID
%         targetRxn  - target reaction ID
%         outputFile - CSV output path (or '' if stdout)
%         matFile    - path to optforce_result.mat
%         rows       - per-intervention rows (grouped by force set)
%         targets    - all must sets + force sets + FVA + diagnostics
    % ------------------- Parameters & inputs -------------------
    if nargin < 6 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end
    % --- Type squash ---
    biomassRxn = char(biomassRxn);
    targetRxn  = char(targetRxn);
    CsourceRxn = char(CsourceRxn);

    % --- Solver setup ---
    % changeCobraSolver('gurobi', 'ALL');
    if nargin < 5
        outputFile = parameters.designDir;
    end
    % --- outputFile default (mirror algFseof/algOptknock) ---
    safeTarget = safeTargetName(targetRxn);
    safeModel  = safeTargetName(modelIdentity(model));
    runBase    = ['OptForce_' safeModel '_' safeTarget];
    csvName    = ['optforce_' safeTarget '.csv'];
    matName    = ['optforce_' safeTarget '_result.mat'];
    outputFile = fullfile(outputFile, csvName);
    algoDir    = fileparts(outputFile);
    matFile    = fullfile(algoDir, matName);
    output = ~isempty(outputFile);

    % --- Re-run cleanup (mirror algFseof/algOptknock) ---
    for oldName = {csvName, matName}
        p = fullfile(algoDir, oldName{1});
        if exist(p, 'file')
            try
                delete(p);
            catch
            end
        end
    end

    % --- Step 1: max growth + max target ---
    model = setParam(model, 'obj', biomassRxn, 1);
    growthRate_sol = solveLP(model);
    fprintf('[algOptforce] The maximum growth rate is %1.2f\n', growthRate_sol.f);

    model = setParam(model, 'obj', targetRxn, 1);
    maxTarget_sol = solveLP(model);
    fprintf('[algOptforce] The maximum production rate of Target is %1.2f\n', maxTarget_sol.f);

    wtBiomass     = growthRate_sol.f;
    targetMaxFlux = maxTarget_sol.f;

    % --- Step 2: WT / mutant constraints ---
    constrWT = struct('rxnList', {{biomassRxn}}, 'rxnValues', 0.98*wtBiomass, 'rxnBoundType', 'b');

    constrMT = struct('rxnList', {{biomassRxn, targetRxn}}, 'rxnValues', [0.01*wtBiomass, 0.98*targetMaxFlux], ...
                      'rxnBoundType', 'bb');

    % --- Step 3: FVA on both strains ---
    [minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ~, ~] = FVAOptForce(model, constrWT, constrMT);

    % --- Step 4a: First-order must sets ---
    runID = [runBase '_K1'];
    constrOpt = struct('rxnList', {{CsourceRxn, biomassRxn, targetRxn}}, 'values', [-10, 0.01*wtBiomass, 0.98*targetMaxFlux]');
    mustLName  = ['MustL_' safeTarget];
    mustUName  = ['MustU_' safeTarget];
    mustUUName = ['MustUU_' safeTarget];
    mustLLName = ['MustLL_' safeTarget];
    mustULName = ['MustUL_' safeTarget];
    forceName  = ['OptForce_' safeTarget];

    [mustLSet, pos_mustL] = findMustL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                                      'runID', runID, 'outputFolder', 'OutputsFindMustL', ...
                                      'outputFileName', mustLName, 'printExcel', 1, 'printText', 1, ...
                                      'printReport', 1, 'keepInputs', 1, 'verbose', 0);

    [mustUSet, pos_mustU] = findMustU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                                      'runID', runID, 'outputFolder', 'OutputsFindMustU', ...
                                      'outputFileName', mustUName, 'printExcel', 1, 'printText', 1, ...
                                      'printReport', 1, 'keepInputs', 1, 'verbose', 0);

    % --- Step 4b: Second-order must sets ---
    exchangeRxns = model.rxns(cellfun(@isempty, strfind(model.rxns, 'EX_')) == 0);
    excludedRxns = unique([mustUSet; mustLSet; exchangeRxns]);

    [mustUU, pos_mustUU, mustUU_linear, pos_mustUU_linear] = ...
        findMustUU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                   'excludedRxns', excludedRxns,'runID', runID, ...
                   'outputFolder', 'OutputsFindMustUU', 'outputFileName', mustUUName, ...
                   'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
                   'verbose', 1);

    [mustLL, pos_mustLL, mustLL_linear, pos_mustLL_linear] = ...
        findMustLL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                   'excludedRxns', excludedRxns,'runID', runID, ...
                   'outputFolder', 'OutputsFindMustLL', 'outputFileName', mustLLName, ...
                   'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
                   'verbose', 1);

    [mustUL, pos_mustUL, mustUL_linear, pos_mustUL_linear] = ...
        findMustUL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustUL', 'outputFileName', mustULName, ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);

    % --- Aggregate must sets ---
    mustU = unique(union(mustUSet, mustUU));
    mustL = unique(union(mustLSet, mustLL));

    % --- Step 5: optForce (K=1, then K=2) ---
    k = 1;
    nSets = 1;
    constrOpt_force = struct('rxnList', {{CsourceRxn, biomassRxn}}, 'values', [-100, 0]);

    [optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
        optForce(model, targetRxn, biomassRxn, mustU, mustL, ...
                 minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
                 'k', k, 'nSets', nSets, 'constrOpt', constrOpt_force, ...
                 'runID', runID, 'outputFolder', 'OutputsOptForce', ...
                 'outputFileName', forceName, 'printExcel', 1, 'printText', 1, ...
                 'printReport', 1, 'keepInputs', 1, 'verbose', 1);

    % --- Second optForce call (K=2, nSets=20) ---
    k = 2;
    nSets = 20;
    runID2 = [runBase '_K2'];
    excludedRxns2 = struct('rxnList', {{optForceSets}}, 'typeReg','U');
    [optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
        optForce(model, targetRxn, biomassRxn, mustU, mustL, ...
                 minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
                 'k', k, 'nSets', nSets, 'constrOpt', constrOpt_force, ...
                 'excludedRxns', excludedRxns2, ...
                 'runID', runID2, 'outputFolder', 'OutputsOptForce', ...
                 'outputFileName', forceName, 'printExcel', 1, 'printText', 1, ...
                 'printReport', 1, 'keepInputs', 1, 'verbose', 1);

    % --- Build result struct (algFseof-style 7 fields) ---
    result = struct( ...
        'config',     struct('K', 2, 'NSets', nSets, ...
                             'WTGrowthFrac', 0.95, 'MTGrowthFrac', 0.5, ...
                             'CSourceBound', -10, ...
                             'CsourceRxn', CsourceRxn, ...
                             'RunID', runID2), ...
        'biomassRxn', biomassRxn, ...
        'targetRxn',  targetRxn, ...
        'outputFile', outputFile, ...
        'matFile',    '', ...
        'rows',       struct('setID', {}, 'interventionType', {}, 'rxnID', {}, ...
                             'rxnName', {}, 'subsystems', {}, 'grRule', {}, ...
                             'postFlux', {}), ...
        'targets',    struct( ...
            'mustUSet',      {{mustUSet}}, ...
            'mustLSet',      {{mustLSet}}, ...
            'mustUU',        {{mustUU}}, ...
            'mustLL',        {{mustLL}}, ...
            'mustUL',        {{mustUL}}, ...
            'mustU',         {{mustU}}, ...
            'mustL',         {{mustL}}, ...
            'forceSets',     {{optForceSets}}, ...
            'typeReg',       {{typeRegOptForceSets}}, ...
            'fluxes',        {flux_optForceSets}, ...
            'minFluxesW',    minFluxesW, ...
            'maxFluxesW',    maxFluxesW, ...
            'minFluxesM',    minFluxesM, ...
            'maxFluxesM',    maxFluxesM, ...
            'wtBiomass',     wtBiomass, ...
            'targetMaxFlux', targetMaxFlux));

    % --- Open CSV (file or stdout) and write header ---
    if output
        outputFile = char(outputFile);
        fid = fopen(outputFile, 'w');
        fprintf(fid, 'SetID,InterventionType,RxnID,RxnName,Subsystems,GrRule,PostFlux\n');
    else
        fprintf('SetID,InterventionType,RxnID,RxnName,Subsystems,GrRule,PostFlux\n');
    end

    % --- Build result.rows and write CSV (per-set, grouped, separated) ---
    nSetsOut = numel(optForceSets);
    for k = 1:nSetsOut
        rxnsInSet   = optForceSets{k};
        typesInSet  = typeRegOptForceSets{k};
        fluxesInSet = flux_optForceSets{k};
        for j = 1:numel(rxnsInSet)
            % --- Build row data ---
            rowFlux = 0;
            if ~isempty(fluxesInSet) && numel(fluxesInSet) >= j
                rowFlux = fluxesInSet(j);
            end
            rxnName    = '';
            subsystems = '';
            grRule     = '';
            pos = find(strcmp(model.rxns, rxnsInSet{j}), 1);
            if ~isempty(pos)
                if isfield(model, 'rxnNames')
                    rxnName = char(model.rxnNames{pos});
                end
                if isfield(model, 'subSystems')
                    ss = model.subSystems{pos};
                    if ~iscell(ss), ss = {ss}; end
                    subsystems = strjoin(ss, ';');
                end
                if isfield(model, 'grRules')
                    grRule = char(model.grRules{pos});
                end
            end

            % --- Append to result.rows ---
            result.rows(end+1) = struct( ...
                'setID',            k, ...
                'interventionType', typesInSet{j}, ...
                'rxnID',            rxnsInSet{j}, ...
                'rxnName',          rxnName, ...
                'subsystems',       subsystems, ...
                'grRule',           grRule, ...
                'postFlux',         rowFlux);

            % --- Write to CSV ---
            if output
                fprintf(fid, '%d,%s,%s,%s,%s,%s,%g\n', ...
                    k, csvEscape(typesInSet{j}), csvEscape(rxnsInSet{j}), ...
                    csvEscape(rxnName), csvEscape(subsystems), csvEscape(grRule), rowFlux);
            else
                fprintf('%d,%s,%s,%s,%s,%s,%g\n', ...
                    k, csvEscape(typesInSet{j}), csvEscape(rxnsInSet{j}), ...
                    csvEscape(rxnName), csvEscape(subsystems), csvEscape(grRule), rowFlux);
            end
        end
        % Separator between sets
        if output
            fprintf(fid, '\n');
        else
            fprintf('\n');
        end
    end

    if output
        fclose(fid);
    end

    % --- Persist result.mat ---
    if ~isempty(algoDir) && ~exist(algoDir, 'dir')
        mkdir(algoDir);
    end
    result.matFile = matFile;
    save(result.matFile, 'result');
    fprintf('[algOptforce] Saved `result` to %s\n', result.matFile);
end


function safe = safeTargetName(targetRxn)
    safe = regexprep(char(targetRxn), '[^A-Za-z0-9_-]+', '_');
    safe = regexprep(safe, '^_+|_+$', '');
    if isempty(safe), safe = 'target'; end
end


function identity = modelIdentity(model)
% modelIdentity  Return a stable model name for OptForce run artifacts.
    for field = {'id', 'name', 'modelID', 'modelName'}
        if isfield(model, field{1}) && ~isempty(model.(field{1}))
            value = model.(field{1});
            if ischar(value) || (isstring(value) && isscalar(value))
                candidate = strtrim(char(value));
                if ~isempty(candidate)
                    identity = candidate;
                    return;
                end
            end
        end
    end
    identity = 'model';
end


function s = csvEscape(s)
% csvEscape  Quote a string for CSV output if it contains a comma, double
%   quote, or newline. Embedded double quotes are doubled per RFC 4180.
    if isempty(s), return; end
    if any(s == ',') || any(s == '"') || any(s == sprintf('\n')) || any(s == sprintf('\r'))
        s = strrep(s, '"', '""');
        s = ['"' s '"'];
    end
end
