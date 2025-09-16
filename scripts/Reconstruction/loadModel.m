function model = loadModel(filename, modeltype, modelDir, parameters)
% LOADMODEL Load and standardize a metabolic model.
%   model = loadModel(filename, modeltype, modelDir, parameters)
%
%   Inputs (all optional):
%     - filename   : filename or path to model file. If empty, parameters.InitialModel is used.
%     - modeltype  : model type string (ECOMAP, TRADITION, SMOMENT, ECMPY, GECKO)
%     - modelDir   : directory to resolve relative filenames (fallback: parameters.modelDir or pwd)
%     - parameters : struct returned by ParameterManager.getParams() (if omitted, it will be fetched)
%
%   The function supports at least .xml, .json, .yml and .yaml extensions. It will
%   attempt to resolve relative filenames against modelDir. Loading and standardizing
%   steps are wrapped to provide clearer error messages.

    % 1) Load parameters once
    if nargin < 4 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.')
        end
    end

    % 2) Resolve modelDir
    if nargin < 3 || isempty(modelDir)
        if isfield(parameters, 'modelDir') && ~isempty(parameters.modelDir)
            modelDir = parameters.modelDir;
        else
            error('The modelDir is not specified');
        end
    end

    % 3) Resolve modeltype
    if nargin < 2 || isempty(modeltype)
        if isfield(parameters, 'modeltype') && ~isempty(parameters.modeltype)
            modeltype = parameters.modeltype;
        else
            error('The modeltype is not specified');
        end
    end
    modeltype = upper(string(modeltype));

    % 4) Resolve filename
    if nargin < 1 || isempty(filename)
        if isfield(parameters, 'InitialModel') && ~isempty(parameters.InitialModel)
            filename = parameters.InitialModel;
        else
            error('loadModel:NoFilename', 'No filename provided and parameters.InitialModel is empty.');
        end
    end

    % If the provided filename is not found as-is, try resolving against modelDir.
    if ~isfile(filename)
        alt = fullfile(modelDir, filename);
        if isfile(alt)
            filename = alt;
        else
            error('loadModel:FileNotFound', 'File not found: %s', filename);
        end
    end

    % 5) Validate modeltype
    validTypes = ["ECOMAP","TRADITION","SMOMENT","ECMPY","GECKO"];
    if ~any(modeltype == validTypes)
        error('loadModel:BadModelType', ...
              'Unsupported model type ''%s''. Valid types: %s', char(modeltype), strjoin(cellstr(validTypes), ', '));
    end

    % 6) Initialize unified model structure
    model = initializeModelStruct();

    % 7) Choose loader based on extension (case-insensitive)
    [~, nameNoExt, ext] = fileparts(filename);
    ext = lower(ext);

    switch ext
        case '.xml'
            loader = @loadXMLmodel;
        case '.json'
            loader = @loadJSONmodel;
        case {'.yml', '.yaml'}
            loader = @loadYMLmodel;
        otherwise
            error('loadModel:BadExt', 'Unsupported file extension: %s', ext);
    end

    % 8) Load file (wrap errors for context)
    try
        model = loader(model, filename);
    catch ME
        % Re-throw with clearer identifier and file context
        error('loadModel:LoadFailed', 'Failed to load model from ''%s'': %s', filename, ME.message);
    end

    % 9) Standardize model structure (wrap errors)
    try
        model = standardizeModel(model, nameNoExt, char(modeltype));
    catch ME
        error('loadModel:StandardizeFailed', 'Failed to standardize model ''%s'': %s', nameNoExt, ME.message);
    end
    model.metNames = regexprep(model.metNames, '\s+(?:\(?[A-Z][a-z]?\d*\)?){2,}(?:[+-]\d*)?$', '');
end