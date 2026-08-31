function result = mdpRunOptforce(model_id, target, biomass, k, nsets, max_candidates)
%MDPRUNOPTFORCE optForce algorithm  --  MUST/U/L and forced interventions.
%   result = mdpRunOptforce(model_id, target, biomass, k, nsets, max_candidates)
%   returns the standard bridge envelope (see CONTRACT.md) carrying:
%       .candidates  --  algorithm-defined result struct (targets + diagnostics)
%       .summary     --  short struct with at least 'count' and the K/NSets
%
%   The bridge defers to ``algOptforce`` (in StrainDesign/algorithms) and
%   catches solver-licence failures with `err_gurobi_license` per the
%   contract catalogue.
    here = fileparts(mfilename('fullpath'));
    addpath_once(here);                                                          % bridge/ for _bridge_util
    addpath_once(fullfile(here, '..', '..', 'StrainDesign'));                    % scripts/StrainDesign
    addpath_once(fullfile(here, '..', '..', 'StrainDesign', 'algorithms'));      % scripts/StrainDesign/algorithms

    if nargin < 1 || isempty(model_id)
        result = make_err('err_param_invalid', 'model_id required');
        return;
    end
    if nargin < 3 || isempty(target) || isempty(biomass)
        if isempty(target),  result = make_err('err_no_target',  'target reaction required'); return; end
        if isempty(biomass), result = make_err('err_no_biomass', 'biomass reaction required'); return; end
    end
    if nargin < 4 || isempty(k),     k     = 2; end
    if nargin < 5 || isempty(nsets), nsets = 1; end
    if nargin < 6 || isempty(max_candidates), max_candidates = 500; end

    try
        model = resolve_model_id(model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end

    bridge_log('mdpRunOptforce', 'target=%s bio=%s K=%d nsets=%d maxCand=%d', ...
               char(target), char(biomass), double(k), double(nsets), double(max_candidates));

    try
        if exist('algOptforce', 'file') ~= 2 && exist('algOptforce', 'file') ~= 5
            error('algOptforce:notFound', 'algOptforce helper not on path');
        end
        candidates = algOptforce(model, char(biomass), char(target), [], '');
    catch err
        result = make_err('err_gurobi_license', err.message);
        return;
    end
    if isempty(candidates) || ~isstruct(candidates)
        candidates = struct();
        count = 0;
    else
        count = 0;
        if isfield(candidates, 'targets') && isstruct(candidates.targets) && numel(fieldnames(candidates.targets)) > 0
            fns = fieldnames(candidates.targets);
            for j = 1:numel(fns)
                v = candidates.targets.(fns{j});
                if ~isempty(v), count = count + numel(v); end
            end
        elseif isfield(candidates, 'rows') && ~isempty(candidates.rows)
            count = numel(candidates.rows);
        end
    end
    payload = struct('candidates', candidates, ...
                     'summary',    struct('count', count, ...
                                         'k', double(k), ...
                                         'nsets', double(nsets), ...
                                         'max_candidates', double(max_candidates)));
    result = make_ok(payload);
end
