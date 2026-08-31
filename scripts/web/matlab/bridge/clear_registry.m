function clear_registry()
%CLEAR_REGISTRY Empty the model registry. Intended for unit tests.
    REG = registry_get();
    remove(REG, keys(REG));   % keep the same Map handle, just drop all entries
end