function model = resolve_model_id(model_id)
%RESOLVE_MODEL_ID Look up a registered model by its UUID-style string id.
% Throws an MException with identifier 'ECOMAP:err_model_format' if the id is
% unknown. Models must be registered via register_model / register_test_model
% before they can be resolved.
    REG = registry_get();
    if ~isKey(REG, model_id)
        throw(MException('ECOMAP:err_model_format', ...
            sprintf('unknown model_id: %s', model_id)));
    end
    model = REG(model_id);
end