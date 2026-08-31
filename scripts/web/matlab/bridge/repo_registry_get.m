function REG = repo_registry_get()
%REPO_REGISTRY_GET Return the shared KcatRepo registry (one persistent containers.Map).
%
%   Mirrors registry_get for the KcatRepo registry used by mdpKcatRepoInit
%   and any future bridge that needs to fetch a previously created repo.
%   See registry_get for the rationale (persistent scoping).
    persistent REPO_REGISTRY
    if isempty(REPO_REGISTRY) || ~isa(REPO_REGISTRY, 'containers.Map')
        REPO_REGISTRY = containers.Map;
    end
    REG = REPO_REGISTRY;
end