function models = BatchModelgeneration(model, P, nutrExch)

    if ~(numel(P) == size(nutrExch,2))
        error('adjBaseModel: incorrect input dimensions')
    end
    
    models = cell(1,numel(P));
    excRxns = nutrExch.Properties.RowNames;
    % set upper and lower bounds for all reactions to zero
    model = setParam(model, 'eq', model.rxns(ismember(model.rxns,excRxns)), 0);
    for i=1:numel(models)
        m = model;
        % set all bounds with available information
        lbVals  = nutrExch{:, i};
        ubVals  = ones(size(lbVals)) * 1000;  % upper bound = 1000
        m = setParam(m, 'lb', excRxns, lbVals);
        m = setParam(m, 'ub', excRxns, ubVals);
        models{i} = m;
    end
end
