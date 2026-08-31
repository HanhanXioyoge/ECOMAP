function REG = registry_get()
%REGISTRY_GET Return the shared model registry (one persistent containers.Map).
%
%   All bridge helper functions that need to store or retrieve a model by
%   id MUST go through this function rather than declaring their own
%   `persistent REGISTRY`. In MATLAB, `persistent` variables are scoped to
%   the single function file that declares them, so each helper used to
%   have its own empty registry and `register_model` would silently drop
%   models that `resolve_model_id` could never see.
%
%   Bridge code outside this folder should call
%   `register_model` / `resolve_model_id` / `clear_registry`, not this.
    persistent REGISTRY
    if isempty(REGISTRY) || ~isa(REGISTRY, 'containers.Map')
        REGISTRY = containers.Map;
    end
    REG = REGISTRY;
end