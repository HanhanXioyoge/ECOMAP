function result = mdpBayesian(ec_model_id, scenarios, max_iter, proc, num_per_gen, reject_num, run_gauks_after)
%MDPBAYESIAN Bayesian MCMC over one or more scenarios.
%   result = mdpBayesian(ec_model_id, scenarios, max_iter, proc, num_per_gen,
%                         reject_num, run_gauks_after) returns the standard
%   bridge envelope (see CONTRACT.md) carrying:
%       .results  --  array of structs, one per scenario, with fields:
%                       .name           scenario tag (e.g. 'ABC_S')
%                       .ec_model_id    UUID of the calibrated ecModel
%                       .rmse_history   per-iteration RMSE
%
%   Per Q2, this function performs all of its work through the single
%   MATLAB Engine session. It does NOT spin up a fresh parpool per call;
%   any parallelism is delegated to the underlying bayesianTuning helper.
    here = fileparts(mfilename('fullpath'));
    addpath_once(here);                                                  % bridge/ for _bridge_util
    addpath_once(fullfile(here, '..', '..', 'Calibration'));             % scripts/Calibration
    addpath_once(fullfile(here, '..', '..', 'Calibration', 'Bayesian')); % scripts/Calibration/Bayesian

    if nargin < 1 || isempty(ec_model_id)
        result = make_err('err_param_invalid', 'ec_model_id required');
        return;
    end
    if nargin < 2 || isempty(scenarios)
        result = make_err('err_param_invalid', 'scenarios list required');
        return;
    end
    if nargin < 3 || isempty(max_iter),       max_iter    = 150; end
    if nargin < 4 || isempty(proc),           proc        = 1;   end
    if nargin < 5 || isempty(num_per_gen),    num_per_gen = 160; end
    if nargin < 6 || isempty(reject_num),     reject_num  = 0.2; end
    if nargin < 7 || isempty(run_gauks_after), run_gauks_after = false; end

    try
        ecModel = resolve_model_id(ec_model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end

    bridge_log('mdpBayesian', 'iter=%d proc=%d nPerGen=%d reject=%.2f postGauks=%d', ...
               double(max_iter), double(proc), double(num_per_gen), ...
               double(reject_num), logical(run_gauks_after));

    out_results = struct('name', {}, 'ec_model_id', {}, 'rmse_history', {});
    useCell = iscell(scenarios);
    nScen = numel(scenarios);
    for k = 1:nScen
        s = scenarios(k);
        if useCell, s = s{1}; end
        if isstruct(s) && numel(s) == 1
            useC = logical(getFieldOr(s, 'use_constraint', false));
            useU = logical(getFieldOr(s, 'use_unconstrained', false));
            useF = logical(getFieldOr(s, 'use_13c_flux', false));
            name = char(getFieldOr(s, 'name', sprintf('scenario_%d', k)));
        else
            % Plain scalar/short mode tag; treat as constraint-only.
            useC = true; useU = false; useF = false;
            name = sprintf('scenario_%d', k);
        end
        bridge_log('mdpBayesian', 'scenario %d/%d "%s" (C=%d U=%d F=%d)', ...
                   k, nScen, name, useC, useU, useF);
        try
            [ecModel_out, ~, ~, rmseHistory, ~] = bayesianTuning( ...
                ecModel, useC, useU, useF, true, ...
                double(max_iter), double(proc), double(num_per_gen), ...
                double(reject_num));
        catch err
            result = make_err('err_gurobi_license', err.message);
            return;
        end
        new_id = char(java.util.UUID.randomUUID.toString);
        register_model(new_id, ecModel_out);
        if run_gauks_after
            bridge_log('mdpBayesian', 'GAUKS follow-up after "%s"', name);
            try
                [ecModel_out2, ~, ~] = GAUKS(ecModel_out, ecModel, 'r_2111');
                new_id = char(java.util.UUID.randomUUID.toString);
                register_model(new_id, ecModel_out2);
            catch err
                bridge_log('mdpBayesian', 'post-GAUKS failed: %s', err.message);
            end
        end
        if isempty(rmseHistory), rmseHistory = []; end
        out_results(k).name         = name;
        out_results(k).ec_model_id  = new_id;
        out_results(k).rmse_history = rmseHistory(:)';
    end
    payload = struct('results', out_results);
    result = make_ok(payload);
end

function v = getFieldOr(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end
