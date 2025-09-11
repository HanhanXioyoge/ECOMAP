function f = objectiveFcn_1(x, ecModel, fit_data, kcat_list)

    % global ecModel;
    % global fit_data;
    % global kcat_list
    reaction_list = kcat_list;
    iterate_model = ecModel;

    objective_value = 0;

    % change kcat
    kcats = x(1:end);
    % kcats = upper_bounds;
    % kcats = start_kcats;

    for i = 1:numel(reaction_list)
        reaction_name = reaction_list{i};
        reaction_index = find(strcmp(reaction_name, iterate_model.ec.rxns));
        new_kcat = kcats(i);
        iterate_model.ec.kcat(reaction_index) = new_kcat;
    end
    
    iterate_model = applyKcatConstraints(iterate_model);

    % Calculate the predicted value
    fba_error = false;
    sum_weights = 0;
    for i = 1:length(fit_data)
        data = fit_data{i}.setup;
        fit_model = fit_data_apply(iterate_model, data);
        sol = solveLP(fit_model);
        predict_value = sol.f;
        if isempty(predict_value)
            disp("FBA error");
            fba_error = true;
        end

        % Calculate objective value
        target_value = fit_data{i}.target;
        objective_value_addition = abs( predict_value / target_value );
        if (objective_value_addition < 1.0)
            objective_value_addition = 1 / objective_value_addition;
        end
        objective_value_addition = objective_value_addition - 1;
        if isfield(fit_data{i}, 'weight')
            objective_value_addition = objective_value_addition * fit_data{i}.weight;
            sum_weights = sum_weights + fit_data{i}.weight;
        else
            sum_weights = sum_weights + 1;
        end
        objective_value = objective_value + abs(objective_value_addition);
    end
    % Return objective value
    objective_value = objective_value / sum_weights;
    if ~fba_error
        f = objective_value;
        disp(f)
    else
        f = 100000;
    end
    % sprintf('%.12f', objective_value)
end

    
