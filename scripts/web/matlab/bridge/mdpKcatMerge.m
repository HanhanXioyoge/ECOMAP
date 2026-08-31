function result = mdpKcatMerge(ec_model_ids, use_custom_file, options, manager_path)
%MDPKCATMERGE Integrate kcat sources and assign values to ecModels.
%   Single-topology runs stay in that topology's reaction naming space.
%   Multi-topology runs build a nameMap and remap kcat lists from the selected
%   reference topology to each target topology only when needed.
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Reconstruction'));
    if nargin < 2 || isempty(use_custom_file)
        use_custom_file = false;
    end
    if nargin < 3 || isempty(options)
        options = struct();
    end
    if nargin >= 4 && ~isempty(manager_path)
        try
            ParameterManager.getParams(manager_path);
        catch err
            result = make_err('err_param_invalid', err.message);
            return;
        end
    end

    ecModel_ids_cell = ec_model_ids;
    if ~iscell(ecModel_ids_cell)
        ecModel_ids_cell = {ecModel_ids_cell};
    end
    if isempty(ecModel_ids_cell)
        result = make_err('err_param_invalid', 'at least one ecModel is required');
        return;
    end

    ecModels = cell(1, numel(ecModel_ids_cell));
    modelTypes = cell(1, numel(ecModel_ids_cell));
    for k = 1:numel(ecModel_ids_cell)
        try
            ecModels{k} = resolve_model_id(ecModel_ids_cell{k});
            modelTypes{k} = ec_model_type(ecModels{k});
        catch err
            result = make_err('err_model_format', err.message);
            return;
        end
    end

    baseType = lower(option_string(options, 'kcatReferenceTopology', choose_base_type(modelTypes)));
    if ~ismember(baseType, modelTypes)
        baseType = choose_base_type(modelTypes);
    end
    baseIdx = find(strcmp(modelTypes, baseType), 1);
    predictionModel = option_string(options, 'kcatPredictionModel', 'CatPred');
    customRxnNameType = lower(option_string(options, 'customKcatRxnNameType', baseType));
    medianThreshold = option_number(options, 'medianThreshold', 5);
    useLoggedMedian = option_bool(options, 'useLoggedMedian', true);

    nameMap = table();
    matchTbl = table();
    mapTbl = table();
    if numel(ecModels) > 1
        bridge_log('mdpKcatMerge', 'Mapping reaction names across %d model type(s)', numel(ecModels));
        try
            [nameMap, matchTbl] = map_ec_rxnname_intersect(ecModels);
            idxBasic = find(strcmp(modelTypes, 'basic'), 1);
            idxIntegrated = find(strcmp(modelTypes, 'integrated'), 1);
            if ~isempty(idxBasic) && ~isempty(idxIntegrated)
                mapTbl = map_basic_one_enzyme_per_row(ecModels{idxBasic}, ecModels{idxIntegrated});
            end
        catch err
            result = make_err('err_kcat_merge', err.message);
            return;
        end
    end

    bridge_log('mdpKcatMerge', 'Retrieving %s kcat predictions from %s model', predictionModel, baseType);
    try
        predictionBase = getPrediction(ecModels{baseIdx}, predictionModel);
        completeBase = completeKcatMatch(ecModels{baseIdx}, struct(), predictionModel);
    catch err
        result = make_err('err_kcat_merge', err.message);
        return;
    end

    mergedIds = struct();
    summary = struct('total', 0, 'from_db', 0, 'from_dl', 0, 'from_median', 0);
    remapStats = struct();
    for i = 1:numel(ecModels)
        mdl = ecModels{i};
        targetType = modelTypes{i};
        pred = predictionBase;
        comp = completeBase;
        if ~strcmpi(targetType, baseType)
            try
                [pred, ~, droppedPred, statsPred] = remapKcatListByType(nameMap, baseType, targetType, predictionBase);
                [comp, ~, droppedComp, statsComp] = remapKcatListByType(nameMap, baseType, targetType, completeBase);
                remapStats.(targetType) = struct( ...
                    'prediction', statsPred, ...
                    'complete', statsComp, ...
                    'predictionDropped', numel(droppedPred), ...
                    'completeDropped', numel(droppedComp));
            catch err
                result = make_err('err_kcat_merge', err.message);
                return;
            end
        end

        try
            fuzzy = fuzzyKcatMatch(mdl);
        catch err
            bridge_log('mdpKcatMerge', 'fuzzyKcatMatch warning: %s', err.message);
            fuzzy = struct();
        end
        try
            merged = mergeKcats(comp, pred, fuzzy);
            mdl = selectKcatValue(mdl, merged);
        catch err
            result = make_err('err_kcat_merge', err.message);
            return;
        end

        if use_custom_file
            try
                [mdl, ~, ~] = fillCustomKcats(mdl, '', customRxnNameType, nameMap);
            catch err
                bridge_log('mdpKcatMerge', 'fillCustomKcats warning: %s', err.message);
            end
        end
        try
            [mdl, ~] = fillMissingKcatWithMedian(mdl, medianThreshold, useLoggedMedian);
        catch err
            bridge_log('mdpKcatMerge', 'fillMissingKcatWithMedian warning: %s', err.message);
        end
        try
            mdl = UpdateSmatrix(mdl);
            mdl = updateProtPool(mdl, true);
        catch err
            bridge_log('mdpKcatMerge', 'updateProtPool warning: %s', err.message);
        end

        newId = char(java.util.UUID.randomUUID.toString);
        register_model(newId, mdl);
        mergedIds.(targetType) = newId;
        if isfield(merged, 'kcat')
            summary.total = summary.total + numel(merged.kcat);
        end
        if isfield(merged, 'source')
            sources = merged.source;
            if ~iscell(sources), sources = num2cell(sources); end
            summary.from_db = summary.from_db + sum(strcmp(sources, 'DB'));
            summary.from_dl = summary.from_dl + sum(strcmp(sources, 'DL'));
            summary.from_median = summary.from_median + sum(strcmp(sources, 'Median'));
        end
    end

    payload = struct( ...
        'merged_ec_model_ids', mergedIds, ...
        'kcat_summary', summary, ...
        'reference_topology', baseType, ...
        'prediction_model', predictionModel, ...
        'remap_stats', remapStats, ...
        'match_count', height(nameMap), ...
        'match_table_count', height(matchTbl), ...
        'basic_integrated_map_count', height(mapTbl));
    result = make_ok(payload);
end

function value = ec_model_type(model)
    value = '';
    if isfield(model, 'enzymeConstraints') && isfield(model.enzymeConstraints, 'ecModeltype')
        value = lower(char(string(model.enzymeConstraints.ecModeltype)));
    end
    if isempty(value)
        error('ecModel is missing enzymeConstraints.ecModeltype.');
    end
end

function value = choose_base_type(modelTypes)
    for candidate = {'integrated', 'isozyme', 'basic'}
        if any(strcmp(modelTypes, candidate{1}))
            value = candidate{1};
            return;
        end
    end
    value = modelTypes{1};
end

function value = option_string(options, field, defaultValue)
    value = defaultValue;
    if isstruct(options) && isfield(options, field) && ~isempty(options.(field))
        value = char(string(options.(field)));
    end
end

function value = option_number(options, field, defaultValue)
    value = defaultValue;
    if isstruct(options) && isfield(options, field) && ~isempty(options.(field))
        parsed = double(options.(field));
        if isfinite(parsed)
            value = parsed;
        end
    end
end

function value = option_bool(options, field, defaultValue)
    value = defaultValue;
    if isstruct(options) && isfield(options, field)
        value = logical(options.(field));
    end
end
