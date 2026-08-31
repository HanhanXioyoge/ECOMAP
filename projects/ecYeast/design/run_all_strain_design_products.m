%% Run selected strain-design algorithms for all 44 target metabolites
% Prerequisites from the model-construction section:
%   GEM
%   ecGEM_integrated
%   ecGEM_isozyme
%   ecGEM_basic
%
% This script replaces the old single-target Step 2 and Step 3 blocks.
% The product loop, directory creation, checkpoints, error isolation and
% summary export are implemented by runEcYeastStrainDesignBatch().
% The batch function uses r_2111 as the biomass reaction and r_1714 as the
% substrate uptake reaction.

%% 1. Validate input models
requiredVariables = { ...
    'GEM', ...
    'ecGEM_integrated', ...
    'ecGEM_isozyme', ...
    'ecGEM_basic'};

missingVariables = requiredVariables(~cellfun( ...
    @(name) evalin('base', ['exist(''' name ''',''var'') == 1']), ...
    requiredVariables));

if ~isempty(missingVariables)
    error('run_all_strain_design_products:MissingModels', ...
        'Load or construct these model variables first: %s', ...
        strjoin(missingVariables, ', '));
end

%% 2. Resolve project parameters and register the batch function
parameters = ParameterManager.getParams();
if isempty(parameters) || ~isstruct(parameters)
    error('run_all_strain_design_products:MissingParameters', ...
        'ParameterManager is not configured for the ecYeast project.');
end

if ~exist('projectDir', 'var') || isempty(projectDir)
    if isfield(parameters, 'projectDir') && ~isempty(parameters.projectDir)
        projectDir = parameters.projectDir;
    else
        projectDir = fullfile(findECOMAProot, 'projects', 'ecYeast');
    end
end
projectDir = char(projectDir);
addpath(fullfile(projectDir, 'design'));

%% 3. Select algorithms
% true  = run this algorithm for every product
% false = do not run this algorithm
options = struct( ...
    ... % Resume behavior
    'SkipCompleted', true, ...
    ...
    ... % GEM algorithms
    'RunFSEOF', true, ...
    'RunOptKnock', true, ...
    'RunOptForce', false, ...
    ...
    ... % ecGEM FSEOF algorithms
    'RunEcFSEOFIntegrated', true, ...
    'RunEcFSEOFIsozyme', true, ...
    'RunEcFSEOFBasic', true, ...
    ...
    ... % OKO algorithms
    'RunOKO', true, ...
    'RunOKOPlus', true, ...
    ...
    ... % COBRA LP/MILP solver
    'CobraSolver', 'gurobi', ...
    ...
    ... % OKO/OKO+ solver settings
    'OkoTimeLimit', 900, ...
    'OkoVerbose', true, ...
    ...
    ... % OKO+ interval source
    'OkoPlusPredictors', {{'unikp'}});

%% 4. Select OKO+ interval sources
% Run UniKP only:
% options.OkoPlusPredictors = {'unikp'};

% Run CatPred only:
% options.OkoPlusPredictors = {'catpred'};

% Run DLKcat only:
% options.OkoPlusPredictors = {'dlkcat'};

% Run all three interval sources:
% options.OkoPlusPredictors = {'unikp', 'catpred', 'dlkcat'};

%% 5. Configure the solver
% The batch function configures the LP solver automatically. It also
% configures the MILP solver when OptKnock or OptForce is enabled.
% Change options.CobraSolver above if a different COBRA solver is required.

%% 6. Run the internal product loop
% Do not add another targetRxn loop around this function call.
% runEcYeastStrainDesignBatch() already iterates over all 44 targets.
batchSummary = runEcYeastStrainDesignBatch( ...
    GEM, ...
    ecGEM_integrated, ...
    ecGEM_isozyme, ...
    ecGEM_basic, ...
    projectDir, ...
    parameters, ...
    options);

%% 7. Show the execution summary
if isempty(batchSummary)
    fprintf('[batch] No algorithm records were produced.\n');
else
    [groupId, algorithm, status] = findgroups( ...
        batchSummary.algorithm, batchSummary.status);
    count = splitapply(@numel, batchSummary.targetRxn, groupId);
    executionSummary = table(algorithm, status, count, ...
        'VariableNames', {'algorithm', 'status', 'productCount'});
    executionSummary = sortrows(executionSummary, {'algorithm', 'status'});
    disp(executionSummary);
end

resultsDir = fullfile(projectDir, 'design', 'strain_design_batch');
summaryFile = fullfile(resultsDir, 'batch_summary.csv');

fprintf('\nResults directory:\n%s\n', resultsDir);
fprintf('\nDetailed summary:\n%s\n', summaryFile);
