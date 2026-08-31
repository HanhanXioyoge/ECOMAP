function result = mdpGrowthPredict(ec_model_id, c_source, bio_rxn, manager_path)
%MDPGROWTHPREDICT Predict growth rate under three conditions; save ecModel.
%   result = mdpGrowthPredict(ec_model_id, c_source, bio_rxn) returns the
%   standard bridge envelope (see CONTRACT.md) carrying:
%       .gem_growth        --  FBA on the original GEM (no enzyme constraints)
%       .unlimited_growth  --  FBA on the ecModel with unlimited protein pool
%       .limited_growth    --  FBA on the ecModel with the limited protein pool
%       .saved_to          --  path to the saved limited ecModel .mat
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
    if nargin < 3 || isempty(bio_rxn)
        result = make_err('err_no_biomass', 'bioRxn not provided');
        return;
    end
    bridge_log('mdpGrowthPredict', 'Solving FBA: original GEM');
    try
        original = setParam(ecModel, 'lb', c_source, -10);
        original = setParam(original, 'obj', bio_rxn, 1);
        sol_gem = solveLP(original);
    catch err
        result = make_err('err_gurobi_license', err.message);
        return;
    end
    bridge_log('mdpGrowthPredict', 'Solving FBA: ecModel with unlimited protein');
    unlimited = setParam(ecModel, 'lb', 'prot_pool_exchange', -1000);
    try
        sol_unlimited = solveLP(unlimited);
    catch err
        result = make_err('err_gurobi_license', err.message);
        return;
    end
    bridge_log('mdpGrowthPredict', 'Solving FBA: ecModel with limited protein');
    try
        limited = updateProtPool(ecModel, true);
        sol_limited = solveLP(limited);
    catch err
        result = make_err('err_gurobi_license', err.message);
        return;
    end
    out_path = fullfile(pwd, [ec_model_id '.mat']);
    try
        save(out_path, 'limited');
    catch err
        bridge_log('mdpGrowthPredict', 'save warning: %s', err.message);
    end
    payload = struct( ...
        'gem_growth',       abs(sol_gem.f), ...
        'unlimited_growth', abs(sol_unlimited.f), ...
        'limited_growth',   abs(sol_limited.f), ...
        'saved_to',         out_path);
    result = make_ok(payload);
end
