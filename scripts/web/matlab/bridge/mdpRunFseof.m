function jsonStr = mdpRunFseof(path, biomassRxn, targetRxn, iterations, coefficient)
%MDPRUNFSEOF Run FSEOF; return {columns, rows, log} in the standard envelope.
%   Failures come back as make_err envelopes (err_no_target / err_no_biomass /
%   err_model_format) rather than raising, so the HTTP layer can report the
%   error code instead of 500ing on a MATLAB exception.
    arguments
        path (1,:) char
        biomassRxn (1,:) char
        targetRxn (1,:) char
        iterations (1,1) double
        coefficient (1,1) double
    end
    log = {};
    try
        model = loadModel(path);
    catch err
        jsonStr = jsonencode(make_err('err_model_format', err.message));
        return;
    end
    if ~ismember(targetRxn, model.rxns)
        jsonStr = jsonencode(make_err('err_no_target', ...
            sprintf('target reaction "%s" not in model', targetRxn)));
        return;
    end
    if ~ismember(biomassRxn, model.rxns)
        jsonStr = jsonencode(make_err('err_no_biomass', ...
            sprintf('biomass reaction "%s" not in model', biomassRxn)));
        return;
    end
    log{end+1} = sprintf('FSEOF start: iters=%d coeff=%.3g', iterations, coefficient);
    try
        r = algFseof(model, biomassRxn, targetRxn, iterations, coefficient);
        tbl = mapFseofResult(r);
    catch err
        jsonStr = jsonencode(make_err('err_param_invalid', err.message));
        return;
    end
    log{end+1} = sprintf('FSEOF done: %d candidates', height(tbl));
    cols = tbl.Properties.VariableNames;
    cell2d = table2cell(tbl);   % MxN cell array
    % jsonencode flattens a 2-D cell array to a 1-D array. Wrap each row in
    % its own cell so the JSON output is a nested list-of-lists.
    nRows = height(tbl);
    rows = cell(1, nRows);
    for i = 1:nRows
        rows{i} = cell2d(i, :);   % 1xN cell (one row)
    end
    out = struct('columns',{cols}, 'rows',{rows}, 'log',{log});
    jsonStr = jsonencode(make_ok(out));
end
