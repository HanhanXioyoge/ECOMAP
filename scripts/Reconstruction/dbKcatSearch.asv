function kcatList = dbKcatSearch(model, parameters)
%%
    % ------- Normalize inputs -------
    if nargin < 2 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.');
        end
    end

    if ~isfield(model,'enzymeConstraints') || ~isfield(model.enzymeConstraints,'eccodes')
        error('No EC codes in model.enzymeConstraints.eccodes.');
    end

    if ~isfield(parameters,'org_name') || isempty(parameters.org_name)
        error('parameters.org_name is required (e.g., ''Escherichia coli'').');
    end

    ecRxns = true(numel(model.enzymeConstraints.rxns),1);
    ecRxns = find(ecRxns);

    % ------- Collect EC codes & substrates for the selected ecRxns -------
    eccodes      = model.enzymeConstraints.eccodes(ecRxns);
    substrates   = cell(numel(ecRxns),1);
    substrCoeffs = cell(numel(ecRxns),1);

    if ~isfield(model.enzymeConstraints, 'ecModeltype') || ~strcmp(model.enzymeConstraints.ecModeltype, 'basic')
        rxnNames = model.enzymeConstraints.rxns;
    else
        rxnNames = extractAfter(model.enzymeConstraints.rxns,4);
    end

    

