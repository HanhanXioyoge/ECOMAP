function result = mdpLoadModel(file_path, model_type, manager_path)
%MDPLOADMODEL Load a metabolic model from disk via the StrainDesign wrapper.
%   result = mdpLoadModel(file_path, model_type, manager_path) returns the standard bridge
%   envelope (see CONTRACT.md) carrying:
%       .model_id           (char)    --  UUID handle for the loaded model
%       .type               (char)    --  'GEM' | 'ecModel'
%       .rxns, .mets, .genes (int)    --  count of each list
%       .detected_organism  (char)    --  value of model.id when present, else ''
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'StrainDesign'));
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Reconstruction'));
    if nargin < 2 || isempty(model_type)
        model_type = 'Tradition';
    end
    parameters = [];
    if nargin >= 3 && ~isempty(manager_path)
        try
            parameters = ParameterManager.getParams(manager_path);
        catch err
            result = make_err('err_model_format', err.message);
            return;
        end
    end
    bridge_log('mdpLoadModel', 'Loading model %s (type=%s)', file_path, model_type);
    try
        model = loadModel(file_path, model_type, [], parameters);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end
    detected = '';
    if isfield(model, 'id') && ~isempty(model.id)
        detected = char(model.id);
    end
    detectedType = 'GEM';
    if isfield(model, 'ecModel') && model.ecModel
        detectedType = 'ecModel';
    end
    new_id = char(java.util.UUID.randomUUID.toString);
    register_model(new_id, model);
    payload = struct( ...
        'model_id', new_id, ...
        'type', detectedType, ...
        'rxns', numel(model.rxns), ...
        'mets', numel(model.mets), ...
        'genes', numel(model.genes), ...
        'detected_organism', detected);
    result = make_ok(payload);
end
