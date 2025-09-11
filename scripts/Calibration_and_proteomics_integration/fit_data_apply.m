function Model = fit_data_apply(model, data)
% fit_data_apply
%   This function is used to characterize the imposition of data environment constraints
%
% Input:
%   model       a model structure
%   data        contains the structures of the reaction that needs to be changed
%               reaction_id:                
%                   'lb': lower bound
%                   'ub': upper bound

% Get the name of the reaction that needs to set
reaction_names = fieldnames(data);
Model = model;
% setting
for i = 1:length(reaction_names)
    reaction_name = reaction_names{i};
    constraints = eval("data."+reaction_name);
    bounds = fieldnames(constraints);
    for j = 1:length(bounds)
        if bounds{j} == "lower_bound"
            lb = constraints.lower_bound;
            Model = setParam(Model, 'lb', reaction_name, lb);
        end
        if bounds{j} == "upper_bound"
            ub = constraints.upper_bound;
            Model = setParam(Model, 'ub', reaction_name, ub);
        end
    end
end
end