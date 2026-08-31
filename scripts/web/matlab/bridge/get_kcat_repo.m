function r = get_kcat_repo(repo_id)
%GET_KCAT_REPO Retrieve a previously registered KcatRepo by id.
    REG = repo_registry_get();
    if ~isKey(REG, repo_id)
        throw(MException('ECOMAP:err_kcat_merge', ...
            sprintf('unknown kcat_repo_id: %s', repo_id)));
    end
    r = REG(repo_id);
end