function result = mdpBuildOkoIntervalsFromModel(model, predictors, options, parameters)
%MDPBUILDOKOINTERVALSFROMMODEL  OKO+ interval pipeline that takes the model directly.
%
%   result = mdpBuildOkoIntervalsFromModel(model, predictors, params)
%   result = mdpBuildOkoIntervalsFromModel(model, predictors, params, options)
%
%   Same semantics as mdpBuildOkoIntervalsFromHomologs but takes the
%   model struct directly — no need to register_model + resolve_model_id.
%   Intended for manual / interactive / scripted use.
%
%   Inputs:
%     model     - ecModel struct (must have .enzymeConstraints and .S).
%     predictors - cellstr from {'UniKP', 'CatPred'}.
%     parameters    - struct from ParameterManager.getParams() OR a char path
%                 to a *ParameterManagement.m file. The function will load
%                 the file if a char is given.
%     options   - optional struct with injection points (same as
%                 runOkoPipeline); default values match the bridge.
%
%   Returns the standard envelope {ok, error_code, error_message, result}
%   where result carries:
%     .predictor_csv_paths      - struct mapping predictor -> CSV path
%     .n_candidates_per_enzyme  - table(rxn, nHomologs)
%     .qc                       - aggregated QC table
%     .elapsed_seconds          - duration
%
%   Example (manual use):
%
%     setup;
%     model = loadModel('eciML1515/models/ecModel_iML1515.mat');
%     parameters = ParameterManager.getParams( ...
%         'eciML1515/eciML1515ParameterManagement.m');
%     result = mdpBuildOkoIntervalsFromModel(model, {'UniKP', 'CatPred'});
%     disp(result.result.predictor_csv_paths);

    % Resolve params struct: accept a path string and load.
    if nargin < 4 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.')
        end
    end

    if nargin < 3 || isempty(options)
        options = struct();
    end

    % Predictors allowlist
    if nargin < 2 || isempty(predictors)
        predictors = {'UniKP', 'CatPred'};
    end
    if ~iscell(predictors)
        predictors = cellstr(predictors);
    end
    allowed = {'UniKP', 'CatPred'};
    bad = setdiff(predictors, allowed);
    if ~isempty(bad)
        result = make_err('err_param_invalid', ...
            sprintf('Unsupported predictors: %s', strjoin(bad, ', ')));
        return;
    end

    % Validate model
    if ~isstruct(model) || ~isfield(model, 'enzymeConstraints') || ...
            ~isfield(model, 'S')
        result = make_err('err_model_format', ...
            'model must be a struct with .enzymeConstraints and .S');
        return;
    end

    % Apply default injection points if not provided
    if ~isfield(options, 'CacheDir') || isempty(options.CacheDir)
        options.CacheDir = fullfile(parameters.reconstructionDir, '.cache', 'uniprot');
    end
    if ~isfield(options, 'dockerFn') || isempty(options.dockerFn)
        options.dockerFn = @defaultDockerFn;
    end
    if ~isfield(options, 'Logger') || isempty(options.Logger)
        options.Logger = @defaultLogger;
    end
    if ~isfield(options, 'organismName') || isempty(options.organismName)
        if isfield(model, 'id') && ~isempty(model.id)
            options.organismName = char(model.id);
        else
            options.organismName = 'ecModel';
        end
        options.organismName = regexprep(options.organismName, '[^\w\-]', '_');
    end
    if ~isfield(options, 'writeDlInput') || isempty(options.writeDlInput)
        options.writeDlInput = @defaultWriteDlInput;
    end
    if ~isfield(options, 'buildSeeds') || isempty(options.buildSeeds)
        options.buildSeeds = @buildOkoReactionSeeds;
    end
    if ~isfield(options, 'fetchEnzymes') || isempty(options.fetchEnzymes)
        options.fetchEnzymes = @fetchCrossSpeciesEnzymes;
    end
    if ~isfield(options, 'buildXlInput') || isempty(options.buildXlInput)
        options.buildXlInput = @buildCrossSpeciesDlInput;
    end
    if ~isfield(options, 'aggregate') || isempty(options.aggregate)
        options.aggregate = @aggregateOkoIntervals;
    end
    if ~isfield(options, 'formatCsv') || isempty(options.formatCsv)
        options.formatCsv = @formatOkoPredictionsCsv;
    end

    % Delegate to the public MATLAB core; this bridge only adds the JSON
    % envelope and never owns OKO+ business logic.
    try
        coreOptions = options;
        coreOptions.Predictors = predictors;
        coreOptions.Parameters = parameters;
        result = buildOkoPlusIntervals(model, coreOptions);
        if ~isfield(result, 'ok')
            % runOkoPipeline returns a plain struct; wrap it.
            result = struct('ok', true, 'error_code', '', ...
                'error_message', '', 'result', result);
        end
    catch err
        result = make_err('err_internal', err.message);
    end
end

% ---------------------------------------------------------------------------
%  Defaults
% ---------------------------------------------------------------------------

function defaultWriteDlInput(theModel, predictor, theParams)
    writeInputFile(theModel, predictor, theParams);
end

function defaultDockerFn(predictor, outputDir)
    ExecutePrediction({predictor}, outputDir);
end

function defaultLogger(varargin)
    prefix = '[mdpBuildOkoIntervalsFromModel] ';
    if numel(varargin) == 1
        fprintf([prefix '%s\n'], varargin{1});
    else
        fmt = ['%s' varargin{1} '\n'];
        args = [{prefix} varargin(2:end)];
        fprintf(fmt, args{:});
    end
end

function r = make_err(code, msg)
    r = struct('ok', false, 'error_code', code, ...
               'error_message', msg, 'result', struct());
end
