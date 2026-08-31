function register_kcat_repo(repo_id, repo)
%REGISTER_KCAT_REPO Store a KcatRepo instance in the in-memory repo registry.
    REG = repo_registry_get();
    REG(repo_id) = repo;
end