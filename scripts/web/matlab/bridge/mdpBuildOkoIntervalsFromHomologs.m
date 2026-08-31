function result = mdpBuildOkoIntervalsFromHomologs(ec_model_id, predictors, parameters, maxHomologs)
%MDPBUILDOKOINTERVALSFROMHOMOLOGS Web bridge: ecModel -> OKO+ interval CSV.
%
%   The first argument accepts EITHER:
%     - char model_id (registered via register_model / register_test_model)
%     - struct model object (with .enzymeConstraints and .S) — bypasses the
%       registry for manual / interactive use.
%
%   If a struct is passed, this function delegates to
%   mdpBuildOkoIntervalsFromModel (same envelope).
%
%   result = mdpBuildOkoIntervalsFromHomologs(ec_model_id, predictors, manager_path)
%
%   Returns the standard bridge envelope (see CONTRACT.md) carrying:
%       .result.predictor_csv_paths.(predictor) -- absolute path to CSV
%       .result.n_candidates_per_enzyme        -- table(rxn, nHomologs)
%       .result.qc                             -- aggregated QC table
%       .result.elapsed_seconds                -- duration in seconds
%
%   Examples:
%
%     % Manual use — model + params, no registry
%     setup;
%     model = loadModel('eciML1515/models/ecModel_iML1515.mat');
%     params = ParameterManager.getParams( ...
%         'eciML1515/eciML1515ParameterManagement.m');
%     r = mdpBuildOkoIntervalsFromHomologs(model, {'UniKP','CatPred'}, params);
%
%     % Web use — register then call with id
%     register_model('eciML1515', model);
%     r = mdpBuildOkoIntervalsFromHomologs('eciML1515', {'UniKP','CatPred'}, '');
    bridgeDir = fileparts(mfilename('fullpath'));
    addpath_once(fullfile(bridgeDir, '..', '..', 'Reconstruction'));
    addpath_once(fullfile(bridgeDir, '..', '..', 'external_kcat_prediction'));

    if nargin < 4 || isempty(maxHomologs), maxHomologs = 100; end

    % If first arg is a struct, delegate to the model-based bridge.
    if isstruct(ec_model_id)
        result = mdpBuildOkoIntervalsFromModel(ec_model_id, predictors, ...
            struct('MaxHomologs', double(maxHomologs)), parameters);
        return;
    end

    if nargin < 3 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.')
        end
    end
    if nargin < 2 || isempty(predictors)
        predictors = {'UniKP', 'CatPred'};
    end
    if ~iscell(predictors)
        predictors = cellstr(predictors);
    end

    % Predictor allowlist
    allowed = {'UniKP', 'CatPred'};
    bad = setdiff(predictors, allowed);
    if ~isempty(bad)
        result = make_err('err_param_invalid', ...
            sprintf('Unsupported predictors: %s', strjoin(bad, ', ')));
        return;
    end

    % Resolve model
    try
        ecModel = resolve_model_id(ec_model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end

    % Resolve parameters (reconstructionDir) for cache + output paths
    if ~isfield(parameters, 'reconstructionDir') || isempty(parameters.reconstructionDir)
        result = make_err('err_param_invalid', ...
            'ParameterManager.reconstructionDir is required');
        return;
    end

    % Bridge-style logger that the runOkoPipeline orchestrator calls.
    bridgeLogger = @(varargin) bridge_pipeline_log(varargin{:});

    try
        runResult = buildOkoPlusIntervals(ecModel, ...
            struct('Predictors', {predictors}, 'Parameters', parameters, ...
                   'Logger', bridgeLogger, 'MaxHomologs', double(maxHomologs)));
    catch err
        errMsg = err.message;
        if contains(errMsg, 'docker', 'IgnoreCase', true)
            result = make_err('err_docker_missing', errMsg);
            return;
        end
        result = make_err('err_kcat_merge', errMsg);
        return;
    end

    payload = struct( ...
        'predictor_csv_paths', runResult.predictor_csv_paths, ...
        'predictor_input_paths', runResult.predictor_input_paths, ...
        'predictor_output_paths', runResult.predictor_output_paths, ...
        'n_candidates_per_enzyme', runResult.n_candidates_per_enzyme, ...
        'qc', runResult.qc, ...
        'run_dir', runResult.run_dir, ...
        'elapsed_seconds', runResult.elapsed_seconds);
    result = make_ok(payload);
end

function bridge_pipeline_log(varargin)
%BRIDGE_PIPELINE_LOG Forward runOkoPipeline log lines to bridge_log.
    if numel(varargin) == 0
        bridge_log('mdpBuildOkoIntervalsFromHomologs', '%s', '');
        return;
    end
    if numel(varargin) == 1
        msg = '';
        if ischar(varargin{1}), msg = varargin{1}; end
        if isstring(varargin{1}), msg = char(varargin{1}); end
        bridge_log('mdpBuildOkoIntervalsFromHomologs', '%s', msg);
        return;
    end
    fmt = varargin{1};
    args = varargin(2:end);
    bridge_log('mdpBuildOkoIntervalsFromHomologs', fmt, args{:});
end
