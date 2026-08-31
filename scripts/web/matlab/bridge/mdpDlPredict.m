function result = mdpDlPredict(ec_model_id, models, manager_path)
%MDPDLPREDICT Generate DL kcat prediction inputs and run predictions (Docker).
%   result = mdpDlPredict(ec_model_id, models) returns the standard bridge
%   envelope (see CONTRACT.md) carrying:
%       .dl_input_paths.(model_name)   --  path to the model's input CSV
%       .dl_result_paths.(model_name)  --  path to the model's output CSV
%   `models` is a cell of DL model names, e.g. {'DLKcat','UniKP','CatPred'}.
%   Note: runtime verification requires a working Docker daemon; the Python
%   TDD path only validates the envelope contract.
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Reconstruction'));
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'external_kcat_prediction'));
    if nargin >= 3 && ~isempty(manager_path)
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
    bridge_log('mdpDlPredict', 'Generating DL inputs for %d model(s)', numel(models));
    try
        Table = writeInputFile(ecModel, models);
    catch err
        result = make_err('err_param_invalid', err.message);
        return;
    end
    bridge_log('mdpDlPredict', 'Executing DL predictions (Docker)');
    try
        ExecutePrediction(models, '');
    catch err
        if ~isempty(strfind(err.message, 'docker')) || ~isempty(strfind(err.message, 'Docker'))
            result = make_err('err_docker_missing', err.message);
            return;
        end
        result = make_err('err_kcat_merge', err.message);
        return;
    end
    payload = struct( ...
        'dl_input_paths', struct(), ...
        'dl_result_paths', struct());
    for k = 1:numel(models)
        m = models{k};
        payload.dl_input_paths.(m)  = fullfile(pwd, sprintf('in_%s.csv',  m));
        payload.dl_result_paths.(m) = fullfile(pwd, sprintf('out_%s.csv', m));
    end
    result = make_ok(payload);
end
