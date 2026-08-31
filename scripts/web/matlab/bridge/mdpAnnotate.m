function result = mdpAnnotate(ec_model_ids, annotation_stages, manager_path, options)
%MDPANNOTATE Annotate ecModels with protein-complex / EC / metabolite info.
%   result = mdpAnnotate(ec_model_ids, annotation_stages) returns the standard
%   bridge envelope (see CONTRACT.md) carrying:
%       .stats  --  {found, applied, noSMILES, noInChIKey, noMNX}
%   `annotation_stages` is a cell of any of {'complex', 'ec', 'met'}.
%   The annotated ecModel is re-registered in the model registry under the
%   same id, so subsequent bridges see the annotated version.
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Reconstruction'));
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'ParameterManagement'));
    if nargin < 4 || isempty(options)
        options = struct();
    end
    if nargin >= 3 && ~isempty(manager_path)
        try
            ParameterManager.getParams(manager_path);
        catch err
            result = make_err('err_param_invalid', err.message);
            return;
        end
    end
    annotation_stages = normalize_stages(annotation_stages);
    runComplexAnnotation = option_bool(options, 'runComplexAnnotation', any(strcmp(annotation_stages, 'complex')));
    runEcAnnotation = option_bool(options, 'runEcAnnotation', any(strcmp(annotation_stages, 'ec')));
    runMetaboliteAnnotation = option_bool(options, 'runMetaboliteAnnotation', true);
    runMetaNetXIntegration = option_bool(options, 'runMetaNetXIntegration', true);
    metaboliteMethods = option_methods(options, 'metaboliteSources', 'ABC');
    bridge_log('mdpAnnotate', 'Annotating %d ecModel(s); stages=%s; complex=%d; ec=%d; met=%d; metSources=%s; mnx=%d', ...
               numel(ec_model_ids), strjoin(annotation_stages, ','), runComplexAnnotation, runEcAnnotation, ...
               runMetaboliteAnnotation, metaboliteMethods, runMetaNetXIntegration);
    stats = struct('found', 0, 'applied', 0, ...
                   'noSMILES', 0, 'noInChIKey', 0, 'noMNX', 0);
    complexInfo = [];
    if runComplexAnnotation
        try
            complexInfo = getComplexdata();
        catch err
            result = make_err('err_raven_notfound', err.message);
            return;
        end
    end
    for k = 1:numel(ec_model_ids)
        id = ec_model_ids{k};
        try
            ecModel = resolve_model_id(id);
        catch err
            result = make_err('err_model_format', err.message);
            return;
        end
        if runComplexAnnotation
            try
                [ecModel, ~, ~, applied] = applyComplexdata(ecModel, complexInfo, false);
                stats.applied = stats.applied + numel(applied);
            catch err
                bridge_log('mdpAnnotate', 'applyComplexdata failed: %s', err.message);
            end
        end
        if runEcAnnotation
            try
                ecModel = getECnumber(ecModel);
            catch err
                bridge_log('mdpAnnotate', 'getECnumber failed: %s', err.message);
            end
        end
        if runMetaboliteAnnotation
            try
                if k == 1
                    [ecModel, noS, noI] = getMetinfo(ecModel, 2, metaboliteMethods);
                else
                    [ecModel, noS, noI] = getMetinfo(ecModel, 2);
                end
                stats.noSMILES = stats.noSMILES + numel(noS);
                stats.noInChIKey = stats.noInChIKey + numel(noI);
            catch err
                bridge_log('mdpAnnotate', 'met annotation failed: %s', err.message);
            end
        end
        if runMetaNetXIntegration
            try
                [ecModel, noM] = addMetMetaNetXID(ecModel, 1);
                stats.noMNX = stats.noMNX + numel(noM);
            catch err
                bridge_log('mdpAnnotate', 'MetaNetX annotation failed: %s', err.message);
            end
        end
        register_model(id, ecModel);
    end
    stats.found = stats.applied;
    result = make_ok(struct('stats', stats));
end

function stages = normalize_stages(stages)
    if ischar(stages) || isstring(stages)
        stages = cellstr(stages);
    end
    for i = 1:numel(stages)
        if strcmp(stages{i}, 'metabolite')
            stages{i} = 'met';
        end
    end
end

function value = option_bool(options, field, defaultValue)
    value = defaultValue;
    if isstruct(options) && isfield(options, field)
        value = logical(options.(field));
    end
end

function methods = option_methods(options, field, defaultValue)
    methods = defaultValue;
    if isstruct(options) && isfield(options, field)
        raw = options.(field);
        if ischar(raw) || isstring(raw)
            methods = upper(char(join(string(raw), '')));
        elseif iscell(raw)
            methods = upper(strjoin(cellfun(@char, raw, 'UniformOutput', false), ''));
        elseif isnumeric(raw) || islogical(raw)
            flags = logical(raw);
            letters = {'A','B','C'};
            methods = strjoin(letters(flags(1:min(numel(flags), 3))), '');
        end
    end
    methods = regexprep(methods, '[^ABC]', '');
end
