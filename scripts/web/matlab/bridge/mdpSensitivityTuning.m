function result = mdpSensitivityTuning(ec_model_id, glc_ex, target_growth, factor, multi)
%MDPSENSITIVITYTUNING Run single-condition (and optionally multi-condition) SA.
%   result = mdpSensitivityTuning(ec_model_id, glc_ex, target_growth, factor, multi)
%   returns the standard bridge envelope (see CONTRACT.md) carrying:
%       .ssa_model_id  --  UUID of the single-condition sensitivity model
%       .tuned_kcats  --  count of tuned kcat values
%       .msa_model_id  --  UUID of the multi-condition sensitivity model ('' if !multi)
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Calibration'));
    try
        ecModel = resolve_model_id(ec_model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end
    bridge_log('mdpSensitivityTuning', 'Single-condition sensitivity analysis');
    try
        ecModel_sSA = setParam(ecModel, 'lb', glc_ex, -1000);
        [ecModel_sSA, tunedKcats] = sensitivityTuning(ecModel_sSA, target_growth, factor);
    catch err
        result = make_err('err_gurobi_license', err.message);
        return;
    end
    ssa_id = char(java.util.UUID.randomUUID.toString);
    register_model(ssa_id, ecModel_sSA);
    msa_id = '';
    if multi
        bridge_log('mdpSensitivityTuning', 'Multi-condition sensitivity analysis');
        try
            ecModel_mSA = MulticonditionsensitivityTuning(ecModel, factor);
            msa_id = char(java.util.UUID.randomUUID.toString);
            register_model(msa_id, ecModel_mSA);
        catch err
            bridge_log('mdpSensitivityTuning', 'multi-condition failed: %s', err.message);
        end
    end
    if isstruct(tunedKcats)
        tunedCount = numel(fieldnames(tunedKcats));
    else
        tunedCount = numel(tunedKcats);
    end
    payload = struct( ...
        'ssa_model_id', ssa_id, ...
        'tuned_kcats', tunedCount, ...
        'msa_model_id', msa_id);
    result = make_ok(payload);
end
