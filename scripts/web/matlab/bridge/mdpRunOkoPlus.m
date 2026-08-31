function result = mdpRunOkoPlus(model_id, biomassRxn, targetRxn, intervalPath, profile)
%MDPRUNOKOPLUS  OKO+ algorithm  --  restrict kcat changes to predicted ranges,
%   fall back to enzyme abundance changes at high cost (paper Methods p.15).
%
%   result = mdpRunOkoPlus(model_id, biomassRxn, targetRxn, intervalPath, profile)
%
%   Inputs:
%     model_id     char -- registered model id (resolved via resolve_model_id)
%     biomassRxn   char -- biomass reaction id
%     targetRxn    char -- target exchange reaction id
%     intervalPath char -- absolute path to a CSV or .mat with kcat intervals
%                          (rxn, uniprot/enzyme, min, max columns). Pass ''
%                          to let the caller supply a table via the wrapper.
%     profile      char -- 'auto' | 'legacy' | 'integrated'
%
%   Returns the standard bridge envelope carrying:
%     .columns / .rows              -- kcat change table
%     .abundanceColumns / .abundanceRows -- enzyme abundance change table
%     .diagnostics                  -- solverStatus, engineeredProduct, ...
%     .config                       -- resolved options (AbundanceWeight, ...)
%     .profile                      -- canonical profile name actually used
%
%   Defers to algOkoPlus in scripts/StrainDesign/algorithms/.
    here = fileparts(mfilename('fullpath'));
    addpath_once(here);
    addpath_once(fullfile(here, '..', '..', 'StrainDesign'));
    addpath_once(fullfile(here, '..', '..', 'StrainDesign', 'algorithms'));

    if nargin < 1 || isempty(model_id)
        result = make_err('err_param_invalid', 'model_id required');
        return;
    end
    if nargin < 3 || isempty(biomassRxn)
        result = make_err('err_no_biomass', 'biomassRxn required');
        return;
    end
    if isempty(targetRxn)
        result = make_err('err_no_target', 'targetRxn required');
        return;
    end
    if nargin < 4 || isempty(intervalPath)
        result = make_err('err_param_invalid', ...
            'OKO+ requires a kcat interval CSV path (rxn, uniprot, min, max).');
        return;
    end
    if nargin < 5 || isempty(profile), profile = 'auto'; end

    try
        model = resolve_model_id(char(model_id));
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end

    bridge_log('mdpRunOkoPlus', 'model=%s target=%s bio=%s intervals=%s profile=%s', ...
               char(model_id), char(targetRxn), char(biomassRxn), ...
               char(intervalPath), char(profile));

    try
        if exist('algOkoPlus', 'file') ~= 2 && exist('algOkoPlus', 'file') ~= 5
            error('algOkoPlus:notFound', 'algOkoPlus helper not on path');
        end
        r = algOkoPlus(model, char(biomassRxn), char(targetRxn), ...
                       char(intervalPath), struct('Profile', char(profile)));
    catch err
        result = classifyOkoError(err);
        return;
    end

    result = make_ok(okoPlusPayload(r));
end

function payload = okoPlusPayload(r)
%OKOPLUSPAYLOAD  Flatten a MATLAB algOkoPlus result into a JSON-friendly payload.
    [cols, rows] = tableToColumns(r.changes);
    payload = struct('columns',   {cols}, ...
                     'rows',      {rows}, ...
                     'diagnostics', r.diagnostics, ...
                     'config',    r.config, ...
                     'profile',   r.profile);
    if isfield(r, 'abundanceChanges') && ~isempty(r.abundanceChanges)
        [acols, arows] = tableToColumns(r.abundanceChanges);
        payload.abundanceColumns = acols;
        payload.abundanceRows    = arows;
    else
        payload.abundanceColumns = {};
        payload.abundanceRows    = {};
    end
end

function [cols, rows] = tableToColumns(T)
    if isempty(T)
        cols = {};
        rows = {};
        return;
    end
    cols = T.Properties.VariableNames;
    cells = table2cell(T);
    n = height(T);
    rows = cell(1, n);
    for i = 1:n, rows{i} = cells(i, :); end
end

function result = classifyOkoError(err)
%CLASSIFYOKOERROR  Map MATLAB error identifiers to bridge-contract codes.
    id = err.identifier;
    msg = err.message;
    if contains(id, 'ReactionNotFound') || contains(msg, 'not found')
        if contains(msg, 'biomass') || contains(msg, 'BIOMASS')
            result = make_err('err_no_biomass', msg);
        else
            result = make_err('err_no_target', msg);
        end
        return;
    end
    if contains(id, 'MissingIntervals') || contains(msg, 'intervals')
        result = make_err('err_kcat_merge', msg);
        return;
    end
    if contains(id, 'Gurobi') || contains(id, 'Solver') || contains(id, 'License')
        result = make_err('err_gurobi_license', msg);
        return;
    end
    if contains(id, 'MissingKcat') || contains(id, 'NoKcatPairs')
        result = make_err('err_model_format', msg);
        return;
    end
    result = make_err('err_param_invalid', msg);
end