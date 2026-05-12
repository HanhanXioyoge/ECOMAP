function [rmse, exp_complete, simulated_out] = rmsecal(ecModel, data, constrain, objective, c_source, rxn2block, org_name)
% RMSECAL: Calculates Net Carbon Flux RMSE with Circuit Breaker logic.
%
% Input:
%   ecModel     - Enzyme-constrained model
%   data        - Table with calibration data (must have ReadRowNames=true when reading)
%   constrain   - true for constrained (fixed uptake), false for unconstrained (max growth)
%   objective   - Biomass reaction name (e.g., 'biomass')
%   c_source    - Carbon source exchange reaction (e.g., 'EX_glc__D_e')
%   rxn2block   - Cell array of reactions to block (optional)
%   org_name    - Organism name for anaerobic model (optional)
%
% Output:
%   rmse        - RMSE value (999 if any experiment failed)
%   exp_complete - Experimental data array
%   simulated_out - Simulated flux values
%
% Example:
%   [rmse, ~, ~] = rmsecal(model_k, growthdata, true, 'biomass', 'EX_glc__D_e', {}, 'ecoli');

    if nargin < 7 || isempty(org_name), org_name = 'ecoli'; end
    if nargin < 6 || isempty(rxn2block), rxn2block = {}; end

    num_exps = height(data);
    vars = data.Properties.VariableNames;

    % Get column indices by name
    sub_col = find(strcmp(vars, 'Substrate'));
    uptake_col = find(strcmp(vars, 'Uptake'));
    biomass_col = find(strcmp(vars, 'biomass'));
    o2_col = find(strcmp(vars, 'EX_o2_e'));
    nh4_col = find(strcmp(vars, 'EX_nh4_e'));
    oxavail_col = find(strcmp(vars, 'OxAvail'));

    % Exchange metabolite columns (ace, eth, gly, pyr, co2, o2, nh4)
    ex_mets = {'EX_ac_e', 'EX_etoh_e', 'EX_glyc_e', 'EX_pyr_e', 'EX_co2_e', 'EX_o2_e', 'EX_nh4_e'};
    ex_cols = zeros(1, length(ex_mets));
    for j = 1:length(ex_mets)
        idx = find(strcmp(vars, ex_mets{j}));
        if ~isempty(idx)
            ex_cols(j) = idx;
        end
    end
    [~, ex_idx] = ismember(ex_mets(ex_cols > 0), ecModel.rxns);
    ex_cols = ex_cols(ex_cols > 0);

    prot_idx = find(strcmp(ecModel.rxns, 'prot_pool_exchange'));

    % Determine if anaerobic
    if ~isempty(oxavail_col)
        is_anaero = strcmp(data{:, oxavail_col}, 'anaerobic') | strcmp(data{:, oxavail_col}, 'limited');
    else
        is_anaero = false(num_exps, 1);
    end

    % Build exp_all: [biomass, ace, eth, glyc, pyr, co2, o2, nh4]
    % Handle missing columns and text 'Nan' values gracefully
    nEx = length(ex_cols);
    exp_all = zeros(num_exps, 1 + nEx);
    exp_all(:, 1) = data{:, biomass_col}; % biomass

    for j = 1:nEx
        col = ex_cols(j);
        if col > 0 && col <= width(data)
            % Get column values and handle text 'Nan'
            colData = data{:, col};
            numericVals = zeros(num_exps, 1);

            for k = 1:num_exps
                v = colData(k);
                if iscell(v), v = v{1}; end

                if ischar(v)
                    vTrim = strtrim(v);
                    if strcmpi(vTrim, 'nan') || isempty(vTrim)
                        numericVals(k) = NaN;
                    else
                        numericVals(k) = str2double(vTrim);
                    end
                elseif isnan(v)
                    numericVals(k) = NaN;
                else
                    numericVals(k) = v;
                end
            end
            exp_all(:, j + 1) = numericVals;
        else
            exp_all(:, j + 1) = NaN;
        end
    end

    simulated_out = zeros(num_exps, 1 + nEx);
    rmse_vec = zeros(num_exps, 1);

    % --- CRITICAL: Initialize fail flag ---
    any_failed = false;

    for i = 1:num_exps
        model_tmp = ecModel;
        sub_name = data{i, sub_col};
        sub_uptake = data{i, uptake_col};
        target_mu = data{i, biomass_col};

        if is_anaero(i)
            model_tmp = anaerobicModel(model_tmp, org_name);
        end

        sub_idx = find(strcmp(model_tmp.rxns, sub_name));
        obj_idx = find(strcmp(model_tmp.rxns, objective));

        % Build idx_complete: [substrate uptake, biomass, ex_mets]
        % For ex_mets, we need the actual model reaction indices
        all_ex_idx = zeros(1, nEx);
        for j = 1:nEx
            metName = ex_mets{j};
            % Try to find exchange reaction
            ex_rxn = ['EX_' strrep(metName(4:end), '_e', '__D_e')]; % Convert to standard exchange format
            idx = find(strcmp(model_tmp.rxns, ex_rxn));
            if isempty(idx)
                idx = find(strcmp(model_tmp.rxns, metName));
            end
            all_ex_idx(j) = idx;
        end
        idx_complete = [sub_idx, obj_idx, all_ex_idx(all_ex_idx > 0)];

        % --- 3. CONSTRAINTS & OBJECTIVES ---
        model_tmp.c(:) = 0;

        % Current data: [biomass, ace, eth, glyc, pyr, co2, o2, nh4]
        current_data = exp_all(i, 2:end); % Skip biomass

        % Apply byproduct/nutrient boundaries
        % ex_mets order: ace, eth, glyc, pyr, co2, o2, nh4
        % Products (ace, eth, glyc, pyr, co2) should have lb = 0
        model_tmp.lb(idx_complete(4:end)) = 0; % Products

        % O2 constraint
        o2_idx_in_all = find(strcmp(ex_mets, 'EX_o2_e'));
        if ~isnan(current_data(o2_idx_in_all))
            model_tmp.lb(all_ex_idx(o2_idx_in_all)) = current_data(o2_idx_in_all);
        else
            model_tmp.lb(all_ex_idx(o2_idx_in_all)) = -1000;
        end

        % NH4 constraint
        nh4_idx_in_all = find(strcmp(ex_mets, 'EX_nh4_e'));
        if ~isnan(current_data(nh4_idx_in_all))
            model_tmp.lb(all_ex_idx(nh4_idx_in_all)) = current_data(nh4_idx_in_all);
        else
            model_tmp.lb(all_ex_idx(nh4_idx_in_all)) = -1000;
        end

        if ~constrain
            %% MODE 1: UNCONSTRAINED (Predict Max Growth)
            model_tmp.lb(strcmp(model_tmp.rxns, c_source)) = 0;
            model_tmp.lb(sub_idx) = -1000;
            model_tmp.c(obj_idx) = 1;
            sol = solveLP(model_tmp);

            % Secondary step: Protein Minimization (pFBA)
            if ~isempty(sol.x) && sol.stat == 1
                model_tmp.lb(obj_idx) = sol.f * 0.99;
                model_tmp.c(:) = 0; model_tmp.c(prot_idx) = 1;
                sol = solveLP(model_tmp);
            end
        else
            %% MODE 2: CONSTRAINED (Fixed Uptake)
            model_tmp.lb(strcmp(model_tmp.rxns, c_source)) = 0;
            model_tmp.lb(sub_idx) = ternary(isnan(sub_uptake), -1000, sub_uptake);
            model_tmp.c(obj_idx) = 1;
            sol = solveLP(model_tmp);

            if ~isempty(sol.x) && sol.stat == 1
                model_tmp.lb(obj_idx) = sol.f * 0.99;
                model_tmp.c(:) = 0; model_tmp.c(prot_idx) = 1;
                sol = solveLP(model_tmp);
            end
        end

        % --- 4. CIRCUIT BREAKER & RMSE CALCULATION ---
        if ~isempty(sol.x) && sol.stat == 1
            valid_mask = ~isnan(exp_all(i, :));
            c_factors = ecModel.excarbon(idx_complete);
            c_factors(c_factors == 0) = 1;

            exp_c_flux = exp_all(i, valid_mask) .* c_factors(valid_mask);
            sim_c_flux = sol.x(idx_complete(valid_mask))' .* c_factors(valid_mask);

            rxnblockidx = ismember(model_tmp.rxns, setdiff(rxn2block, model_tmp.rxns(idx_complete(2))));
            sim_block_c = sol.x(rxnblockidx)' .* ecModel.excarbon(rxnblockidx);

            % Unified Net Carbon Flux RMSE
            rmse_vec(i) = sqrt(mean(([exp_c_flux, 0] - [sim_c_flux, sum(sim_block_c)]).^2));
            simulated_out(i, :) = sol.x(idx_complete)';
        else
            % INFEASIBLE: Set flag and break the experiment loop immediately
            any_failed = true;
            break;
        end
    end

    % --- 5. FINAL AGGREGATION ---
    if any_failed
        rmse = 999; % Penalty for total failure
        simulated_out(:, :) = NaN;
    else
        rmse = mean(rmse_vec);
    end
    exp_complete = exp_all;
end

function y = ternary(cond, a, b), if cond, y = a; else, y = b; end, end
