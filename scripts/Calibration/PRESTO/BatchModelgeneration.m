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
        m = setParam(m, 'lb', excRxns, lbVals);
        models{i} = m;
    end
end
