function register_model(model_id, model)
%REGISTER_MODEL Store a model in the in-memory registry under a stable id.
    REG = registry_get();
    REG(model_id) = model;
end