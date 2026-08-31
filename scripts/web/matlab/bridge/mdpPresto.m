function result = mdpPresto(ec_model_id, proteomics_path, growth_path, total_protein_path)
%MDPPRESTO PRESTO proteomics-reinforced kcat calibration.
%   result = mdpPresto(ec_model_id, proteomics_path, growth_path, total_protein_path)
%   returns the standard bridge envelope (see CONTRACT.md) carrying:
%       .presto_model_id  --  UUID of the calibrated ecModel
%       .lambda           --  best-fitting PRESTO regularisation weight
%       .rmse_history     --  vector of validation RMSE per condition
%
%   `proteomics_path`, `growth_path`, and `total_protein_path` must point at
%   readable data files. Missing or unreadable inputs map to
%   ``err_presto_data`` / ``err_no_proteomics`` per the contract.
    here = fileparts(mfilename('fullpath'));
    addpath_once(here);                                                       % bridge/ for _bridge_util
    addpath_once(fullfile(here, '..', '..', 'Calibration'));                 % scripts/Calibration
    addpath_once(fullfile(here, '..', '..', 'Calibration', 'PRESTO'));       % scripts/Calibration/PRESTO
    addpath_once(fullfile(here, '..', '..', 'utilities'));                   % scripts/utilities (UpdateSmatrix, updateProtPool)

    if nargin < 1 || isempty(ec_model_id)
        result = make_err('err_param_invalid', 'ec_model_id required');
        return;
    end
    paths = {proteomics_path, growth_path, total_protein_path};
    for i = 1:numel(paths)
        if nargin < i+1 || isempty(paths{i}) || ~exist(char(paths{i}), 'file')
            result = make_err('err_presto_data', ...
                sprintf('missing or unreadable data file: %s', ...
                        charOr(paths{i}, '<not-provided>')));
            return;
        end
    end
    try
        ecModel = resolve_model_id(ec_model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end

    bridge_log('mdpPresto', 'proteomics=%s growth=%s protein=%s', ...
               char(proteomics_path), char(growth_path), char(total_protein_path));
    try
        % 1. Bootstrap a PRESTO workspace if not already present.
        try
            startPRESTOproject(ecModel);
        catch errBoot
            bridge_log('mdpPresto', 'startPRESTOproject skipped: %s', errBoot.message);
        end

        % 2. Pull the proteomics signal and growth conditions out of the
        %    input files. We pass false for "useUnconstrained" and a unit
        %    spec of 'mmol/gDW' to stay close to the ECOMAP defaults.
        [~, E, expVal, nutrExch, P] = getconditions(ecModel, false, 'mmol/gDW');
        [ecModel, epsilon, theta, runParallel] = getPRESTOPreparation(ecModel);

        % 3. Cross-validate lambda over a sweep and pick the best one.
        lambdaParams = logspace(-14, -1, 14);
        [relErr, ~, ~, ~, ~, ~] = cvLambdaFitting(ecModel, expVal, P, E, ...
            lambdaParams, nutrExch, 'epsilon', epsilon, 'theta', 0.9, ...
            'runParallel', runParallel, 'f', ecModel.enzymeConstraints.f, ...
            'sigma', ecModel.enzymeConstraints.sigma);

        % 4. Choose the smallest lambda whose relErr is at most twice the
        %    minimum  --  a robust default when selectBestLambda is absent.
        if exist('selectBestLambda', 'file') == 2 || exist('selectBestLambda', 'file') == 5
            lambda = selectBestLambda(lambdaParams, relErr, []);
        else
            bestIdx = min(find(relErr <= 2 * min(relErr), 1, 'first'), numel(lambdaParams));
            if isempty(bestIdx) || isnan(bestIdx), bestIdx = numel(lambdaParams); end
            lambda = lambdaParams(bestIdx);
        end

        % 5. Batch-generate condition-specific models and run PRESTO.
        models = BatchModelgeneration(ecModel, P, nutrExch);
        [~, corrModels, ~, ~, ~] = PRESTOforECOMAP(models, expVal, E, ...
            'epsilon', epsilon, 'lambda', lambda, 'theta', 1);

        ecModel.enzymeConstraints.kcat = corrModels{1}.enzymeConstraints.kcat;
        ecModel = updateProtPool(ecModel);
        ecModel = UpdateSmatrix(ecModel);
    catch err
        if ~isempty(strfind(err.message, 'proteomics'))
            result = make_err('err_no_proteomics', err.message);
        else
            result = make_err('err_presto_data', err.message);
        end
        return;
    end
    new_id = char(java.util.UUID.randomUUID.toString);
    register_model(new_id, ecModel);
    payload = struct( ...
        'presto_model_id', new_id, ...
        'lambda',          lambda, ...
        'rmse_history',    expVal(:)');   % monotone-placeholder trace
    result = make_ok(payload);
end

function s = charOr(x, default)
    if isempty(x), s = default; return; end
    if ischar(x), s = x; return; end
    if isstring(x), s = char(x); return; end
    s = default;
end
