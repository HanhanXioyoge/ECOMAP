function result = mdpConvertecModel(model_id, topology, manager_path)
%MDPCONVERTECMODEL Expand a GEM into one or more ecModel topologies.
%   result = mdpConvertecModel(model_id, topology) returns the standard
%   bridge envelope (see CONTRACT.md) carrying:
%       .ec_model_ids.(basic|isozyme|integrated)  --  UUID for the converted model
%       .stats.(basic|isozyme|integrated)         --  {rxns, mets, genes}
%   `topology` is one of 'basic', 'isozyme', 'integrated', 'all'.
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Reconstruction'));
    if nargin >= 3 && ~isempty(manager_path)
        try
            ParameterManager.getParams(manager_path);
        catch err
            result = make_err('err_param_invalid', err.message);
            return;
        end
    end
    if ~any(strcmp(topology, {'basic', 'isozyme', 'integrated', 'all'}))
        result = make_err('err_param_invalid', ['unknown topology: ' topology]);
        return;
    end
    if strcmp(topology, 'all')
        types = {'basic', 'isozyme', 'integrated'};
    else
        types = {topology};
    end
    try
        model = resolve_model_id(model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end
    bridge_log('mdpConvertecModel', 'Converting model %s to ecModel (requested topologies: %s)', ...
               model_id, strjoin(types, ','));
    ec_model_ids = struct();
    stats = struct();
    for i = 1:numel(types)
        t = types{i};
        bridge_log('mdpConvertecModel', '  topology = %s', t);
        try
            ecModel = convertecModel(model, t);
        catch err
            result = make_err('err_model_format', err.message);
            return;
        end
        new_id = char(java.util.UUID.randomUUID.toString);
        register_model(new_id, ecModel);
        ec_model_ids.(t) = new_id;
        stats.(t) = struct('rxns', numel(ecModel.rxns), ...
                            'mets', numel(ecModel.mets), ...
                            'genes', numel(ecModel.genes));
    end
    payload = struct('ec_model_ids', ec_model_ids, 'stats', stats);
    result = make_ok(payload);
end
