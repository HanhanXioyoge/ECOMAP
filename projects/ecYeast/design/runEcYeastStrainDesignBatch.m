function summary = runEcYeastStrainDesignBatch( ...
        GEM, ecGEM_integrated, ecGEM_isozyme, ecGEM_basic, ...
        projectDir, parameters, options)
%RUNECYEASTSTRAINDESIGNBATCH Run all strain-design methods for all products.
%
% Example:
%   addpath(fullfile(projectDir, 'design'));
%   options = struct( ...
%       'SkipCompleted', true, ...
%       'RunFSEOF', true, ...
%       'RunOptKnock', true, ...
%       'RunOptForce', true, ...
%       'RunEcFSEOFIntegrated', true, ...
%       'RunEcFSEOFIsozyme', true, ...
%       'RunEcFSEOFBasic', true, ...
%       'RunOKO', true, ...
%       'RunOKOPlus', true, ...
%       'CobraSolver', 'gurobi', ...
%       'OkoTimeLimit', 900, ...
%       'OkoVerbose', true, ...
%       'OkoPlusPredictors', {{'unikp'}});
%   summary = runEcYeastStrainDesignBatch( ...
%       GEM, ecGEM_integrated, ecGEM_isozyme, ecGEM_basic, ...
%       projectDir, ParameterManager.getParams(), options);
%
% Output layout:
%   <projectDir>/design/strain_design_batch/
%     batch_summary.csv
%     01_r_1549_R_R_-2_3-butanediol_stoichiometric/
%       target_info.csv
%       GEM/FSEOF/
%       GEM/OptKnock/
%       GEM/OptForce/
%       ecGEM/ecFSEOF_integrated/
%       ecGEM/ecFSEOF_isozyme/
%       ecGEM/ecFSEOF_basic/
%       ecGEM/OKO/
%       ecGEM/OKOPlus/UniKP/

    if nargin < 5 || isempty(projectDir)
        error('runEcYeastStrainDesignBatch:MissingProjectDir', ...
            'projectDir is required.');
    end
    if nargin < 6 || isempty(parameters)
        parameters = ParameterManager.getParams();
    end
    if ~isstruct(parameters)
        error('runEcYeastStrainDesignBatch:MissingParameters', ...
            'A valid ParameterManager parameter struct is required.');
    end
    if nargin < 7 || isempty(options)
        options = struct();
    end
    options = applyDefaults(options);
    configureCobraSolvers(options.CobraSolver, ...
        options.RunOptKnock || options.RunOptForce);

    biomassRxn = 'r_2111';
    substrateRxn = 'r_1714';
    parameters.bioRxn = biomassRxn;

    targets = targetTable();
    outputRoot = fullfile(projectDir, 'design', 'strain_design_batch');
    ensureDirectory(outputRoot);
    writetable(targets, fullfile(outputRoot, 'targets.csv'));

    if options.RunFSEOF || options.RunOptKnock || options.RunOptForce
        validateReaction(GEM, biomassRxn, 'GEM biomass');
    end
    if options.RunOptForce
        validateReaction(GEM, substrateRxn, 'GEM carbon source');
    end
    if options.RunOptKnock
        selectedRxns = geneAssociatedReactions(GEM);
    else
        selectedRxns = GEM.rxns;
    end

    if options.RunEcFSEOFIntegrated || options.RunOKO || options.RunOKOPlus
        fprintf('[batch] Preparing integrated ecGEM...\n');
        ecGEM_integrated = prepareEcModel( ...
            ecGEM_integrated, 'integrated', substrateRxn);
    end
    if options.RunEcFSEOFIsozyme
        fprintf('[batch] Preparing isozyme ecGEM...\n');
        ecGEM_isozyme = prepareEcModel( ...
            ecGEM_isozyme, 'isozyme', substrateRxn);
    end
    if options.RunEcFSEOFBasic
        fprintf('[batch] Preparing basic ecGEM...\n');
        ecGEM_basic = prepareEcModel(ecGEM_basic, 'basic', substrateRxn);
    end

    records = cell(0, 8);
    for targetIndex = 1:height(targets)
        targetRxn = targets.targetRxn{targetIndex};
        metabolite = targets.metabolite{targetIndex};
        classification = targets.ecomapClassification{targetIndex};
        productFolder = sprintf('%02d_%s_%s_%s', ...
            targetIndex, targetRxn, safeName(metabolite), classification);
        productDir = fullfile(outputRoot, productFolder);
        ensureDirectory(productDir);
        writetable(targets(targetIndex, :), ...
            fullfile(productDir, 'target_info.csv'));

        fprintf('\n[batch] [%d/%d] %s (%s, %s)\n', ...
            targetIndex, height(targets), metabolite, targetRxn, classification);

        % ---------------- GEM algorithms ----------------
        stageDir = fullfile(productDir, 'GEM', 'FSEOF');
        [status, message, elapsed] = runConfiguredStage( ...
            options.RunFSEOF, hasReaction(GEM, targetRxn), ...
            'Target reaction is missing from GEM.', ...
            stageDir, options.SkipCompleted, ...
            @() algFseof(GEM, biomassRxn, targetRxn, 21, 0.99, ...
                stageDir, parameters));
        records(end+1, :) = resultRecord(targetRxn, metabolite, classification, ...
            'FSEOF', status, elapsed, stageDir, message); %#ok<AGROW>
        writeSummary(records, outputRoot);

        stageDir = fullfile(productDir, 'GEM', 'OptKnock');
        [status, message, elapsed] = runConfiguredStage( ...
            options.RunOptKnock, hasReaction(GEM, targetRxn), ...
            'Target reaction is missing from GEM.', ...
            stageDir, options.SkipCompleted, ...
            @() algOptknock(GEM, biomassRxn, targetRxn, selectedRxns, ...
                1, 15, 0.40, 'G', true, stageDir, parameters));
        records(end+1, :) = resultRecord(targetRxn, metabolite, classification, ...
            'OptKnock', status, elapsed, stageDir, message); %#ok<AGROW>
        writeSummary(records, outputRoot);

        stageDir = fullfile(productDir, 'GEM', 'OptForce');
        [status, message, elapsed] = runConfiguredStage( ...
            options.RunOptForce, hasReaction(GEM, targetRxn), ...
            'Target reaction is missing from GEM.', ...
            stageDir, options.SkipCompleted, ...
            @() runOptForceInDirectory(GEM, biomassRxn, targetRxn, ...
                substrateRxn, stageDir, parameters));
        records(end+1, :) = resultRecord(targetRxn, metabolite, classification, ...
            'OptForce', status, elapsed, stageDir, message); %#ok<AGROW>
        writeSummary(records, outputRoot);

        % ---------------- ecGEM algorithms ----------------
        ecStages = { ...
            'ecFSEOF_integrated', ecGEM_integrated, options.RunEcFSEOFIntegrated; ...
            'ecFSEOF_isozyme', ecGEM_isozyme, options.RunEcFSEOFIsozyme; ...
            'ecFSEOF_basic', ecGEM_basic, options.RunEcFSEOFBasic};
        for stageIndex = 1:size(ecStages, 1)
            algorithm = ecStages{stageIndex, 1};
            ecModel = ecStages{stageIndex, 2};
            enabled = ecStages{stageIndex, 3};
            stageDir = fullfile(productDir, 'ecGEM', algorithm);
            [status, message, elapsed] = runConfiguredStage( ...
                enabled, hasReaction(ecModel, targetRxn), ...
                'Target reaction is missing from this ecGEM variant.', ...
                stageDir, options.SkipCompleted, ...
                @() algecFSEOF(ecModel, targetRxn, substrateRxn, ...
                    16, true, stageDir, parameters));
            records(end+1, :) = resultRecord(targetRxn, metabolite, classification, ...
                algorithm, status, elapsed, stageDir, message); %#ok<AGROW>
            writeSummary(records, outputRoot);
        end

        stageDir = fullfile(productDir, 'ecGEM', 'OKO');
        okoOptions = struct( ...
            'Profile', 'integrated', ...
            'TimeLimit', options.OkoTimeLimit, ...
            'Verbose', options.OkoVerbose, ...
            'OutputFile', fullfile(stageDir, 'oko_kcat_changes.csv'), ...
            'MatFile', fullfile(stageDir, 'oko_result.mat'));
        [status, message, elapsed] = runConfiguredStage( ...
            options.RunOKO, hasReaction(ecGEM_integrated, targetRxn), ...
            'Target reaction is missing from integrated ecGEM.', ...
            stageDir, options.SkipCompleted, ...
            @() algOko(ecGEM_integrated, biomassRxn, targetRxn, okoOptions));
        records(end+1, :) = resultRecord(targetRxn, metabolite, classification, ...
            'OKO', status, elapsed, stageDir, message); %#ok<AGROW>
        writeSummary(records, outputRoot);

        for predictorIndex = 1:numel(options.OkoPlusPredictors)
            predictor = lower(char(options.OkoPlusPredictors{predictorIndex}));
            predictorLabel = predictorDisplayName(predictor);
            intervalCsv = fullfile(projectDir, 'reconstruction', 'kcatData', ...
                'oko_plus_runs', 'ecYeast9_h25_dlkcat_unikp_catpred', ...
                ['ecYeast9_' predictor '_kcat_preds.csv']);
            stageDir = fullfile(productDir, 'ecGEM', 'OKOPlus', predictorLabel);
            plusOptions = struct( ...
                'Profile', 'integrated', ...
                'TimeLimit', options.OkoTimeLimit, ...
                'Verbose', options.OkoVerbose, ...
                'OutputFile', fullfile(stageDir, ...
                    'okoplus_kcat_changes.csv'), ...
                'MatFile', fullfile(stageDir, 'okoplus_result.mat'));
            available = hasReaction(ecGEM_integrated, targetRxn) && isfile(intervalCsv);
            unavailableMessage = okoPlusUnavailableMessage( ...
                ecGEM_integrated, targetRxn, intervalCsv);
            [status, message, elapsed] = runConfiguredStage( ...
                options.RunOKOPlus, available, unavailableMessage, ...
                stageDir, options.SkipCompleted, ...
                @() algOkoPlus(ecGEM_integrated, biomassRxn, ...
                    targetRxn, intervalCsv, plusOptions));
            algorithm = ['OKOPlus_' predictorLabel];
            records(end+1, :) = resultRecord(targetRxn, metabolite, classification, ...
                algorithm, status, elapsed, stageDir, message); %#ok<AGROW>
            writeSummary(records, outputRoot);
        end
    end

    summary = writeSummary(records, outputRoot);
    save(fullfile(outputRoot, 'batch_summary.mat'), 'summary', 'targets', 'options');
    fprintf('\n[batch] Finished. Summary: %s\n', ...
        fullfile(outputRoot, 'batch_summary.csv'));
end


function options = applyDefaults(options)
    defaults = struct( ...
        'SkipCompleted', true, ...
        'RunFSEOF', true, ...
        'RunOptKnock', true, ...
        'RunOptForce', true, ...
        'RunEcFSEOFIntegrated', true, ...
        'RunEcFSEOFIsozyme', true, ...
        'RunEcFSEOFBasic', true, ...
        'RunOKO', true, ...
        'RunOKOPlus', true, ...
        'CobraSolver', 'gurobi', ...
        'OkoTimeLimit', 900, ...
        'OkoVerbose', true, ...
        'OkoPlusPredictors', {{'unikp'}});
    names = fieldnames(defaults);
    for index = 1:numel(names)
        name = names{index};
        if ~isfield(options, name) || isempty(options.(name))
            options.(name) = defaults.(name);
        end
    end
    logicalFields = { ...
        'SkipCompleted', 'RunFSEOF', 'RunOptKnock', 'RunOptForce', ...
        'RunEcFSEOFIntegrated', 'RunEcFSEOFIsozyme', ...
        'RunEcFSEOFBasic', 'RunOKO', 'RunOKOPlus', 'OkoVerbose'};
    for index = 1:numel(logicalFields)
        name = logicalFields{index};
        validateattributes(options.(name), {'logical', 'numeric'}, ...
            {'scalar'}, mfilename, ['options.' name]);
        options.(name) = logical(options.(name));
    end
    validateattributes(options.OkoTimeLimit, {'numeric'}, ...
        {'scalar', 'positive', 'finite'}, mfilename, 'options.OkoTimeLimit');
    if ischar(options.OkoPlusPredictors) || isstring(options.OkoPlusPredictors)
        options.OkoPlusPredictors = cellstr(options.OkoPlusPredictors);
    end
    if ~iscell(options.OkoPlusPredictors)
        error('runEcYeastStrainDesignBatch:Predictors', ...
            'options.OkoPlusPredictors must be a cell array of predictor names.');
    end
    if isstring(options.CobraSolver)
        if ~isscalar(options.CobraSolver)
            error('runEcYeastStrainDesignBatch:CobraSolver', ...
                'options.CobraSolver must be a scalar string or character vector.');
        end
        options.CobraSolver = char(options.CobraSolver);
    end
    if ~ischar(options.CobraSolver) || isempty(strtrim(options.CobraSolver)) || ...
            size(options.CobraSolver, 1) ~= 1
        error('runEcYeastStrainDesignBatch:CobraSolver', ...
            'options.CobraSolver must be a nonempty character vector.');
    end
end


function configureCobraSolvers(solverName, requireMILP)
    [lpOK, lpInstalled] = changeCobraSolver(solverName, 'LP', 1);
    if ~isscalar(lpOK) || ~lpOK
        error('runEcYeastStrainDesignBatch:LPSolver', ...
            ['Could not configure COBRA LP solver "%s" ' ...
             '(installed=%s).'], solverName, mat2str(lpInstalled));
    end

    if requireMILP
        [milpOK, milpInstalled] = changeCobraSolver( ...
            solverName, 'MILP', 1);
        if ~isscalar(milpOK) || ~milpOK
            error('runEcYeastStrainDesignBatch:MILPSolver', ...
                ['Could not configure COBRA MILP solver "%s" ' ...
                 '(installed=%s). OptKnock and OptForce require a ' ...
                 'working MILP solver.'], ...
                solverName, mat2str(milpInstalled));
        end
    end

    lpSolver = getCobraSolver('LP', false);
    validateSolverValue(lpSolver, 'CBT_LP_SOLVER');
    if requireMILP
        milpSolver = getCobraSolver('MILP', false);
        validateSolverValue(milpSolver, 'CBT_MILP_SOLVER');
    end
end


function validateSolverValue(value, variableName)
    if ~ischar(value) || isempty(value) || size(value, 1) ~= 1
        error('runEcYeastStrainDesignBatch:InvalidSolverValue', ...
            ['COBRA global %s must contain one solver name as a ' ...
             'character vector. Run changeCobraSolver before solving.'], ...
            variableName);
    end
end


function targets = targetTable()
    rows = { ...
        'r_1549', '(R,R)-2,3-butanediol', 'native', 'alcohol', 'stoichiometric', 'stoichiometric'; ...
        'r_1683', 'choline', 'native', 'alkaloid', 'kinetic', 'kinetic'; ...
        'r_1916', 'laurate', 'native', 'fatty acids and lipids', 'kinetic', 'kinetic'; ...
        'r_2189', 'oleate', 'native', 'fatty acids and lipids', 'stoichiometric', 'stoichiometric'; ...
        'r_2188', 'hexanoate', 'native', 'organic acid', 'kinetic', 'kinetic'; ...
        'r_2033', 'pyruvate', 'native', 'organic acid', 'stoichiometric', 'stoichiometric'};
    targets = cell2table(rows, ...
        'VariableNames', {'targetRxn', 'metabolite', 'origin', ...
        'chemicalClass', 'ecFactoryClassification', ...
        'ecomapClassification'});
end


function selectedRxns = geneAssociatedReactions(model)
    if ~isfield(model, 'grRules') || numel(model.grRules) ~= numel(model.rxns)
        error('runEcYeastStrainDesignBatch:MissingGeneRules', ...
            'GEM.grRules must align with GEM.rxns.');
    end
    rules = strtrim(string(model.grRules(:)));
    valid = strlength(rules) > 0 & rules ~= "()" & rules ~= "[]";
    selectedRxns = model.rxns(valid);
end


function model = prepareEcModel(model, profile, substrateRxn)
    [model, ~, ~] = fillCustomKcats(model, '', profile);
    model = UpdateSmatrix(model);
    validateReaction(model, substrateRxn, [profile ' carbon source']);
    model.lb(strcmp(model.rxns, substrateRxn)) = -10;
end


function [status, message, elapsed] = runConfiguredStage( ...
        enabled, available, unavailableMessage, stageDir, skipCompleted, callback)
    if ~logical(enabled)
        status = 'disabled';
        message = '';
        elapsed = 0;
        fprintf('[batch]   disabled: %s\n', stageDir);
        return;
    end
    if ~available
        ensureDirectory(stageDir);
        status = 'failed';
        message = unavailableMessage;
        elapsed = 0;
        errorInfo = struct( ...
            'identifier', 'runEcYeastStrainDesignBatch:Unavailable', ...
            'message', unavailableMessage, ...
            'report', unavailableMessage, ...
            'startedAt', char(datetime('now', ...
                'Format', 'yyyy-MM-dd HH:mm:ss')), ...
            'elapsedSeconds', 0);
        save(fullfile(stageDir, 'batch_error.mat'), 'errorInfo');
        fprintf(2, '[batch]   unavailable: %s\n[batch]   %s\n', ...
            stageDir, unavailableMessage);
        return;
    end
    [status, message, elapsed] = runStage(stageDir, skipCompleted, callback);
end


function [status, message, elapsed] = runStage(stageDir, skipCompleted, callback)
    ensureDirectory(stageDir);
    completedFile = fullfile(stageDir, 'batch_completed.mat');
    errorFile = fullfile(stageDir, 'batch_error.mat');
    if skipCompleted && isfile(completedFile)
        status = 'skipped_completed';
        message = '';
        elapsed = 0;
        fprintf('[batch]   skip completed: %s\n', stageDir);
        return;
    end

    startedAt = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    timerValue = tic;
    try
        callback();
        elapsed = toc(timerValue);
        checkpoint = struct('completed', true, 'startedAt', startedAt, ...
            'elapsedSeconds', elapsed);
        save(completedFile, 'checkpoint');
        if isfile(errorFile), delete(errorFile); end
        status = 'completed';
        message = '';
        fprintf('[batch]   completed in %.1f s: %s\n', elapsed, stageDir);
    catch exception
        elapsed = toc(timerValue);
        status = 'failed';
        message = exception.message;
        errorInfo = struct( ...
            'identifier', exception.identifier, ...
            'message', exception.message, ...
            'report', getReport(exception, 'extended', 'hyperlinks', 'off'), ...
            'startedAt', startedAt, ...
            'elapsedSeconds', elapsed);
        save(errorFile, 'errorInfo');
        fprintf(2, '[batch]   failed: %s\n[batch]   %s\n', ...
            stageDir, exception.message);
    end
end


function message = okoPlusUnavailableMessage(model, targetRxn, intervalCsv)
    if ~hasReaction(model, targetRxn)
        message = 'Target reaction is missing from integrated ecGEM.';
    elseif ~isfile(intervalCsv)
        message = ['Interval CSV not found: ' intervalCsv];
    else
        message = '';
    end
end


function result = runOptForceInDirectory( ...
        model, biomassRxn, targetRxn, substrateRxn, stageDir, parameters)
    oldDirectory = pwd;
    cleanup = onCleanup(@() cd(oldDirectory));
    cd(stageDir);
    result = algOptforce(model, biomassRxn, targetRxn, ...
        substrateRxn, stageDir, parameters);
end


function row = resultRecord(targetRxn, metabolite, ecomapClassification, ...
        algorithm, status, elapsed, outputDir, message)
    row = {targetRxn, metabolite, ecomapClassification, algorithm, status, ...
        elapsed, outputDir, message};
end


function summary = writeSummary(records, outputRoot)
    summary = cell2table(records, 'VariableNames', { ...
        'targetRxn', 'metabolite', 'ecomapClassification', 'algorithm', ...
        'status', 'elapsedSeconds', 'outputDir', 'message'});
    writetable(summary, fullfile(outputRoot, 'batch_summary.csv'));
end


function validateReaction(model, reactionId, role)
    if ~hasReaction(model, reactionId)
        error('runEcYeastStrainDesignBatch:ReactionNotFound', ...
            '%s reaction "%s" is missing.', role, reactionId);
    end
end


function tf = hasReaction(model, reactionId)
    tf = isfield(model, 'rxns') && any(strcmp(model.rxns, reactionId));
end


function ensureDirectory(path)
    if ~exist(path, 'dir')
        [ok, message] = mkdir(path);
        if ~ok
            error('runEcYeastStrainDesignBatch:CreateDirectory', ...
                'Could not create "%s": %s', path, message);
        end
    end
end


function safe = safeName(value)
    safe = regexprep(char(value), '[^A-Za-z0-9_-]+', '_');
    safe = regexprep(safe, '^_+|_+$', '');
    if isempty(safe), safe = 'product'; end
end


function name = predictorDisplayName(predictor)
    switch lower(predictor)
        case 'unikp'
            name = 'UniKP';
        case 'catpred'
            name = 'CatPred';
        case 'dlkcat'
            name = 'DLKcat';
        otherwise
            error('runEcYeastStrainDesignBatch:Predictor', ...
                'Unsupported OKO+ predictor "%s".', predictor);
    end
end
