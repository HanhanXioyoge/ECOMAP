function startPRESTOproject(model, parameters)

    if nargin < 2 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end

    if ~isfield(model,'enzymeConstraints') || isempty(model.enzymeConstraints)
        error('The model does not contain the enzymeConstraints structure.');
    end
    
    path = parameters.path;
    templateFile = fullfile(findECOMAProot, 'scripts', 'Calibration', 'PRESTO', 'PRESTOtemplate.m');
    templateText = fileread(templateFile);

    repl.org_name        = parameters.org_name;
    repl.id              = model.id;
    repl.f               = model.enzymeConstraints.f;
    repl.sigma           = model.enzymeConstraints.sigma;
    repl.runParallel     = parameters.PRESTO.runParallel;
    repl.ncpu            = parameters.PRESTO.ncpu;
    repl.nIter           = parameters.PRESTO.nIter;
    repl.epsilon         = parameters.PRESTO.epsilon;
    repl.lambda          = parameters.PRESTO.lambda;
    repl.theta           = parameters.PRESTO.theta;
    
    if parameters.PRESTO.runParallel
        repl.runParallel = 'true';
    else
        repl.runParallel = 'false';
    end

    fields = fieldnames(repl);
    for i = 1:numel(fields)
        key         = fields{i};
        placeholder = ['${' key '}'];
        value       = repl.(key);
        if ~ischar(value)
            value = num2str(value);
        end
        templateText = strrep(templateText, placeholder, value);
    end
    outFile = fullfile(path, [model.id 'PRESTOConfiguration.m']);
    fid = fopen(outFile, 'w');
    fwrite(fid, templateText);
    fclose(fid);
end