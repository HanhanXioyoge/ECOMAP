function result = mdpRunOko(model_id, biomassRxn, targetRxn, profile)
%MDPRUNOKO  OKO algorithm  --  minimise kcat modifications at fixed enzyme allocation.
%
%   result = mdpRunOko(model_id, biomassRxn, targetRxn, profile)
%
%   Inputs:
%     model_id   char  -- registered model id (resolved via resolve_model_id)
%     biomassRxn char  -- biomass reaction id
%     targetRxn  char  -- target exchange reaction id
%     profile    char  -- 'auto' | 'legacy' | 'integrated'
%
%   Returns the standard bridge envelope (see CONTRACT.md) carrying:
%     .columns / .rows        -- kcat change table
%     .diagnostics            -- solverStatus, engineeredProduct, ...
%     .config                 -- resolved options (AbundanceWeight, etc.)
%     .profile                -- canonical profile name actually used
%
%   Defers to algOko in scripts/StrainDesign/algorithms/.
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
    if nargin < 4 || isempty(profile), profile = 'auto'; end

    try
        model = resolve_model_id(char(model_id));
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end

    bridge_log('mdpRunOko', 'model=%s target=%s bio=%s profile=%s', ...
               char(model_id), char(targetRxn), char(biomassRxn), char(profile));

    try
        if exist('algOko', 'file') ~= 2 && exist('algOko', 'file') ~= 5
            error('algOko:notFound', 'algOko helper not on path');
        end
        r = algOko(model, char(biomassRxn), char(targetRxn), ...
                   struct('Profile', char(profile)));
    catch err
        result = classifyOkoError(err);
        return;
    end

    result = make_ok(okoPayload(r));
end

function payload = okoPayload(r)
%OKOPAYLOAD  Flatten a MATLAB algOko result into the JSON-friendly payload.
    T = r.changes;
    [cols, rows] = tableToColumns(T);
    payload = struct('columns',   {cols}, ...
                     'rows',      {rows}, ...
                     'diagnostics', r.diagnostics, ...
                     'config',    r.config, ...
                     'profile',   r.profile);
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