function model = loadModel(filename, modeltype, modelDir, parameters)
% LOADMODEL Load and standardize a metabolic model.
%   model = loadModel(filename, modeltype, modelDir, parameters)
%
%   Inputs (all optional):
%     - filename   : filename or path to model file. If empty, parameters.InitialModel is used.
%     - modeltype  : model type string (ECOMAP, TRADITION, SMOMENT, ECMPY, GECKO)
%     - modelDir   : directory to resolve relative filenames (fallback: parameters.modelsDir or pwd)
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
        if isfield(parameters, 'modelsDir') && ~isempty(parameters.modelsDir)
            modelDir = parameters.modelsDir;
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
    model = loader(model, filename);
    %{
    try
        model = loader(model, filename);
    catch ME
        % Re-throw with clearer identifier and file context
        error('loadModel:LoadFailed', 'Failed to load model from ''%s'': %s', filename, ME.message);
    end
    %}
    % 9) Standardize model structure (wrap errors)
    try
        model.id = nameNoExt;
        model.type = char(modeltype);
    catch ME
        error('loadModel:StandardizeFailed', 'Failed to standardize model ''%s'': %s', nameNoExt, ME.message);
    end
    
    % ---- Ensure Miriam cross-reference containers exist and match lengths ----
    model.metMiriams  = ensureMiriamCell(model, 'metMiriams',  numel(model.mets));
    model.rxnMiriams  = ensureMiriamCell(model, 'rxnMiriams',  numel(model.rxns));
    model.geneMiriams = ensureMiriamCell(model, 'geneMiriams', numel(model.genes));

    model.metNames = regexprep(model.metNames, '\s+(?:\(?[A-Z][a-z]?\d*\)?){2,}(?:[+-]\d*)?$', '');

    % ---- Standardize metabolite names for GECKO models ----
    % Remove compartment suffix like [c], [e], [p] from metabolite names
    model.mets = standardizeMetNames(model.mets, char(modeltype));

    % ---- Standardize reaction names for ECMPY and GECKO models ----
    model.rxns = standardizeRxnNames(model.rxns, char(modeltype));
end

function metNames = standardizeMetNames(metNames, modeltype)
% STANDARDIZEMETNAMES Remove compartment suffix like [c], [e], [p] from metabolite names.

    if strcmp(modeltype, 'GECKO')
        % Remove [x] suffix where x is a single letter (c, e, p, etc.)
        metNames = regexprep(metNames, '\[[a-z]\]$', '');
    end
end

function rxnNames = standardizeRxnNames(rxnNames, modeltype)
% STANDARDIZERNXNNAMES Standardize reaction names for ECMPY and GECKO models.

    if strcmp(modeltype, 'ECMPY')
        % Handle multiple patterns:
        % _reverse_num2 -> _REV_EXP_2
        % _num2 -> _EXP_2
        % _reverse -> _REV
        rxnNames = regexprep(rxnNames, '_reverse_num(\d+)$', '_REV_EXP_$1');
        rxnNames = regexprep(rxnNames, '_num(\d+)$', '_EXP_$1');
        rxnNames = regexprep(rxnNames, '_reverse$', '_REV');

    elseif strcmp(modeltype, 'GECKO')
        % No2 -> _EXP_2, No3 -> _EXP_3 (all reactions ending with No, not _REV)
        hasNo = ~cellfun(@isempty, regexp(rxnNames, 'No\d+$'));
        hasREV = endsWith(rxnNames, '_REV');

        idx = hasNo & ~hasREV;
        rxnNames(idx) = regexprep(rxnNames(idx), 'No(\d+)$', '_EXP_$1');
    end

    % ---- Cleanup orphan _EXP_N suffixes ----
    % If A_EXP_1 exists but A_EXP_2, A_EXP_3, etc. do not exist, remove _EXP_1
    rxnNames = cleanupOrphanExpSuffix(rxnNames);
end

function rxnNames = cleanupOrphanExpSuffix(rxnNames)
% Remove _EXP_1 when it's the only suffix for that base reaction.

    tokenData = regexp(rxnNames, '^(.+)_EXP_(\d+)$', 'tokens');
    expIdx = find(~cellfun(@isempty, tokenData));

    if isempty(expIdx)
        return;
    end

    baseCounts = containers.Map();
    for i = expIdx'
        baseName = tokenData{i}{1}{1};
        if isKey(baseCounts, baseName)
            baseCounts(baseName) = baseCounts(baseName) + 1;
        else
            baseCounts(baseName) = 1;
        end
    end

    for i = expIdx'
        baseName = tokenData{i}{1}{1};
        expNum = str2double(tokenData{i}{1}{2});
        if expNum == 1 && baseCounts(baseName) == 1
            rxnNames{i} = baseName;
        end
    end
end

function C = ensureMiriamCell(model, fieldName, N)
% ensureMiriamCell
%   Guarantees model.(fieldName) is an N×1 cell array.
%   Each cell holds a struct array with fields .name / .value (or empty).
%
% Behavior:
%   - If the field is missing or empty -> initialize to N×1 with empty structs.
%   - If present, must be a cell array with exactly N elements (column-shaped).
%   - Normalizes empties to an empty struct array (struct('name',{},'value',{})).

    emptyEntry = struct('name',{},'value',{});

    if N == 0
        % No items to host; return empty cell to avoid length mismatch errors.
        C = cell(0,1);
        return;
    end

    if ~isfield(model, fieldName) || isempty(model.(fieldName))
        C = repmat({emptyEntry}, N, 1);
        return;
    end

    C = model.(fieldName);
    if ~iscell(C)
        error('%s must be a cell array (got %s).', fieldName, class(C));
    end
    if numel(C) ~= N
        error('%s length (%d) must match target length (%d).', fieldName, numel(C), N);
    end
    C = C(:); % force column

    % Normalize element types: allow empty or struct; coerce empty to emptyEntry
    for i = 1:N
        if isempty(C{i})
            C{i} = emptyEntry;
        elseif ~isstruct(C{i})
            error('%s{%d} must be a struct (or empty).', fieldName, i);
        end
    end
end