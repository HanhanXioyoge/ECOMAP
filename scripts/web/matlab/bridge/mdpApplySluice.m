function result = mdpApplySluice(ec_model_id, ex_rxn_list)
%MDPAPPLYSLUICE Apply the Sluice structure (prerequisite for calibration).
%   result = mdpApplySluice(ec_model_id, ex_rxn_list) returns the standard
%   bridge envelope (see CONTRACT.md) carrying:
%       .sluice_config          --  struct returned by applySluiceStructure
%       .modified_ec_model_id   --  UUID handle for the modified ecModel
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Calibration'));
    try
        ecModel = resolve_model_id(ec_model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end
    if isempty(ex_rxn_list)
        result = make_err('err_sluice_data', 'ex_rxn_list is empty');
        return;
    end
    bridge_log('mdpApplySluice', 'Applying Sluice structure over %d exchange reaction(s)', ...
               numel(ex_rxn_list));
    try
        [ecModel, sluiceConfig] = applySluiceStructure(ecModel, ex_rxn_list);
    catch err
        result = make_err('err_sluice_data', err.message);
        return;
    end
    new_id = char(java.util.UUID.randomUUID.toString);
    register_model(new_id, ecModel);
    payload = struct('sluice_config', sluiceConfig, 'modified_ec_model_id', new_id);
    result = make_ok(payload);
end
