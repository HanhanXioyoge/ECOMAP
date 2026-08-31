function result = mdpKcatCompare(ec_model_id, dl_models, complex_names, manager_path)
%MDPKCATCOMPARE Compare DL-predicted kcat values against database kcat values.
%   result = mdpKcatCompare(ec_model_id, dl_models, complex_names) returns
%   the standard bridge envelope (see CONTRACT.md) carrying:
%       .eval_table   --  struct array {model, rmse, r2}
%       .match_path   --  path to the saved MATCH struct
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'DLmode_evaluation'));
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Reconstruction'));
    if nargin >= 4 && ~isempty(manager_path)
        try
            ParameterManager.getParams(manager_path);
        catch err
            result = make_err('err_param_invalid', err.message);
            return;
        end
    end
    try
        ecModel = resolve_model_id(ec_model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end
    bridge_log('mdpKcatCompare', 'Building kcat matches');
    try
        MATCH = BuildKcatMatches(dl_models, complex_names);
    catch err
        result = make_err('err_kcat_merge', err.message);
        return;
    end
    bridge_log('mdpKcatCompare', 'Analyzing kcat matches');
    try
        OUT = AnalyzeKcatMatches(MATCH, dl_models, true, true);
    catch err
        OUT = struct();
    end
    save_path = fullfile(pwd, sprintf('%s_match.mat', ec_model_id));
    try
        save(save_path, 'MATCH');
    catch err
        result = make_err('err_kcat_merge', err.message);
        return;
    end
    eval_table = struct('model', {}, 'rmse', {}, 'r2', {});
    if isstruct(OUT) && numel(fieldnames(OUT)) > 0
        for k = 1:numel(dl_models)
            eval_table(k).model = dl_models{k};
            if isfield(OUT, dl_models{k})
                eval_table(k).rmse = OUT.(dl_models{k}).rmse;
                eval_table(k).r2   = OUT.(dl_models{k}).r2;
            end
        end
    end
    payload = struct('eval_table', eval_table, 'match_path', save_path);
    result = make_ok(payload);
end
