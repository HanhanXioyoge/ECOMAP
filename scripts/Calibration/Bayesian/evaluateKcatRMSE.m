function [rmse_final, exp_out, sim_out] = evaluateKcatRMSE(ecModel, growthdata, UnconstrainedGrowth, C13Fluxdata, bioRxn, c_source, rxn2block, org_name, C13ReactionMap)
% EVALUATEKCATRSE Calculate RMSE using calibrated kcat values in ecModel
%   Directly uses ecModel.enzymeConstraints.kcat without UpdateSmatrix
%   Carbon sources are closed once at the beginning

    exp_out = []; sim_out = [];

    % Find all carbon-containing exchange reactions with active uptake and close them
    [EXrxn, EXrxnIdx] = getExchangeRxns(ecModel);
    carbonSourcesToClose = [];
    for i = 1:numel(EXrxnIdx)
        rxnIdx = EXrxnIdx(i);
        metIdx = find(ecModel.S(:, rxnIdx));
        if ~isempty(metIdx)
            metFormula = ecModel.metFormulas{metIdx(1)};
            % Carbon is C followed by a number (e.g., CO2, CH4, C6H12O6)
            % Exclude Ca, Cu, Co, Cl, etc. where C is followed by a letter
            if ~isempty(regexp(metFormula, 'C\d', 'once')) && ecModel.lb(rxnIdx) <= -1000
                carbonSourcesToClose = [carbonSourcesToClose, rxnIdx];
            end
        end
    end
    % Close all carbon sources at once
    ecModel.lb(carbonSourcesToClose) = 0;

    % --- Scenario 1: Constrained Growth Experiments ---
    rmse_1 = NaN; exp_1 = []; sim_1 = [];
    if ~isempty(growthdata)
        [rmse_1, exp_1, sim_1] = rmsecal(ecModel, growthdata, true, bioRxn, c_source, rxn2block, org_name);
    end

    % --- Scenario 2: Unconstrained Growth Experiments ---
    rmse_2 = NaN; exp_2 = []; sim_2 = [];
    if ~isempty(UnconstrainedGrowth)
        [rmse_2, exp_2, sim_2] = rmsecal(ecModel, UnconstrainedGrowth, false, bioRxn, c_source, rxn2block, org_name);
    end

    % --- Scenario 3: 13C Flux Data ---
    rmse_3 = NaN;
    if ~isempty(C13Fluxdata) && isstruct(C13Fluxdata) && isfield(C13Fluxdata, 'reactions')
        try
            nConditions = length(C13Fluxdata.conditions);
            rmse_conditions = zeros(nConditions, 1);
            for condIdx = 1:nConditions
                [rmse_c, ~, ~] = rmsecal_13C(ecModel, C13Fluxdata, condIdx, bioRxn, c_source, C13ReactionMap);
                rmse_conditions(condIdx) = rmse_c;
            end
            rmse_3 = mean(rmse_conditions(~isnan(rmse_conditions)));
        catch ME
            warning('13C RMSE calculation failed: %s', ME.message);
            rmse_3 = NaN;
        end
    end

    % --- Weighted Average RMSE ---
    n1 = height(growthdata);
    n2 = height(UnconstrainedGrowth);
    n3 = 0;
    if ~isempty(C13Fluxdata) && isstruct(C13Fluxdata) && isfield(C13Fluxdata, 'reactions')
        n3 = length(C13Fluxdata.reactions);
    end

    rmse_values = []; weights = [];
    if ~isnan(rmse_1) && n1 > 0, rmse_values = [rmse_values, rmse_1]; weights = [weights, n1]; end
    if ~isnan(rmse_2) && n2 > 0, rmse_values = [rmse_values, rmse_2]; weights = [weights, n2]; end
    if ~isnan(rmse_3) && n3 > 0, rmse_values = [rmse_values, rmse_3]; weights = [weights, n3]; end

    if isempty(rmse_values)
        rmse_final = NaN;
    else
        rmse_final = sum(weights .* rmse_values) / sum(weights);
    end

    exp_out = [exp_1; exp_2];
    sim_out = [sim_1; sim_2];
end


function [rmse, exp_complete, simulated_out] = rmsecal(ecModel, data, constrain, objective, c_source, rxn2block, org_name)
% RMSECAL: Calculates Net Carbon Flux RMSE with Circuit Breaker logic.
%   Carbon sources are already closed by the caller (evaluateKcatRMSE)

    num_exps = height(data);
    constraints = data.Properties.VariableNames;
    ex_mets = constraints(4:10); % ace, eth, gly, pyr, co2, o2, nh4
    [~, ex_idx] = ismember(ex_mets, ecModel.rxns);
    prot_idx = find(strcmp(ecModel.rxns, 'prot_pool_exchange'));

    exp_all = table2array(data(:, 2:10));
    is_anaero = strcmp(data{:, 11}, 'anaerobic') | strcmp(data{:, 11}, 'limited');

    simulated_out = zeros(num_exps, 9);
    rmse_vec = zeros(num_exps, 1);
    any_failed = false;

    for i = 1:num_exps
        model_tmp = ecModel;
        sub_name = data{i, 1};
        sub_uptake = data{i, 2};
        target_mu = data{i, 3};

        if is_anaero(i), model_tmp = anaerobicModel(model_tmp, org_name); end

        sub_idx = find(strcmp(model_tmp.rxns, sub_name));
        obj_idx = find(strcmp(model_tmp.rxns, objective));
        idx_complete = [sub_idx, obj_idx, ex_idx];

        model_tmp.c(:) = 0;
        idx_products = ex_idx(1:5);
        idx_nutrients = ex_idx(6:7);
        current_data = exp_all(i, 3:9);

        model_tmp.lb(idx_products) = 0;
        if ~isnan(current_data(6))
            model_tmp.lb(idx_nutrients(1)) = current_data(6);
        else
            model_tmp.lb(idx_nutrients(1)) = -1000;
        end
        if ~isnan(current_data(7))
            model_tmp.lb(idx_nutrients(2)) = current_data(7);
        else
            model_tmp.lb(idx_nutrients(2)) = -1000;
        end

        if ~constrain
            % MODE 1: UNCONSTRAINED (Predict Max Growth)
            % Carbon sources already closed, only open specified substrate
            model_tmp.lb(sub_idx) = -1000;
            model_tmp.c(obj_idx) = 1;
            sol = solveLP(model_tmp);

            if ~isempty(sol.x) && sol.stat == 1
                model_tmp.lb(obj_idx) = sol.f * 0.99;
                model_tmp.c(:) = 0; 
                % model_tmp.c(:) = -1;
                model_tmp.c(prot_idx) = 1;
                sol = solveLP(model_tmp);
            end
        else
            % MODE 2: CONSTRAINED (Fixed Uptake)
            % Carbon sources already closed, only open specified substrate
            model_tmp.lb(sub_idx) = ternary(isnan(sub_uptake), -1000, sub_uptake);
            model_tmp.c(obj_idx) = 1;
            sol = solveLP(model_tmp);

            if ~isempty(sol.x) && sol.stat == 1
                model_tmp.lb(obj_idx) = sol.f * 0.99;
                model_tmp.c(:) = 0; 
                model_tmp.c(:) = -1;
                % model_tmp.c(prot_idx) = 1;
                sol = solveLP(model_tmp);
            end
        end

        if ~isempty(sol.x) && sol.stat == 1
            valid_mask = ~isnan(exp_all(i, :));
            c_factors = ecModel.excarbon(idx_complete);
            c_factors(c_factors == 0) = 1;

            exp_c_flux = exp_all(i, valid_mask) .* c_factors(valid_mask);
            sim_c_flux = sol.x(idx_complete(valid_mask))' .* c_factors(valid_mask);

            rxnblockidx = ismember(model_tmp.rxns, setdiff(rxn2block, model_tmp.rxns(idx_complete(2))));
            sim_block_c = sol.x(rxnblockidx)' .* ecModel.excarbon(rxnblockidx);

            rmse_vec(i) = sqrt(mean(([exp_c_flux, 0] - [sim_c_flux, sum(sim_block_c)]).^2));
            simulated_out(i, :) = sol.x(idx_complete)';
        else
            any_failed = true;
            break;
        end
    end

    if any_failed
        rmse = 999;
        simulated_out(:, :) = NaN;
    else
        rmse = mean(rmse_vec);
    end
    exp_complete = exp_all;
end

function y = ternary(cond, a, b), if cond, y = a; else, y = b; end, end