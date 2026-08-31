function result = mdpKcatRepoInit(ec_model_id)
%MDPKCATREPOINIT Initialise a fresh KcatRepo with the model's 'Init' kcat set.
%   result = mdpKcatRepoInit(ec_model_id) returns the standard bridge
%   envelope (see CONTRACT.md) carrying:
%       .kcat_repo_id  --  UUID handle for the new KcatRepo
%       .init_group    --  'Init'
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Calibration'));
    try
        ecModel = resolve_model_id(ec_model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end
    if isfield(ecModel, 'information') && isfield(ecModel.information, 'organism')
        organism = ecModel.information.organism;
    else
        organism = 'unknown';
    end
    if isfield(ecModel, 'id') && ~isempty(ecModel.id)
        modelName = ecModel.id;
    else
        modelName = ec_model_id;
    end
    try
        kcatRepo = KcatRepo(modelName, organism);
    catch err
        result = make_err('err_kcat_merge', err.message);
        return;
    end
    if isfield(ecModel, 'enzymeConstraints') && isfield(ecModel.enzymeConstraints, 'kcat')
        Init_kcats = ecModel.enzymeConstraints.kcat;
    else
        Init_kcats = [];
    end
    try
        ex_rxn_list = evalin('base', 'ex_rxn_list');
    catch
        ex_rxn_list = {};
    end
    try
        sluiceParams = extractSluiceParams(ecModel, ex_rxn_list);
    catch err
        sluiceParams = struct();
    end
    try
        kcatRepo.addGroup('Init', Init_kcats, sluiceParams, ...
                          'Deep learning + databases', 'initial kcat');
    catch err
        result = make_err('err_kcat_merge', err.message);
        return;
    end
    repo_id = char(java.util.UUID.randomUUID.toString);
    register_kcat_repo(repo_id, kcatRepo);
    payload = struct('kcat_repo_id', repo_id, 'init_group', 'Init');
    result = make_ok(payload);
end
