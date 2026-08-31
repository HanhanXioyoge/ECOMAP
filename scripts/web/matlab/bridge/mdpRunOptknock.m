function result = mdpRunOptknock(model_id, target, biomass, max_candidates, num_del, min_growth_fraction, vmax)
%MDPRUNOPTKNOCK OptKnock algorithm  --  enumerate reaction-knockout sets.
%   result = mdpRunOptknock(model_id, target, biomass, max_candidates,
%                           num_del, min_growth_fraction, vmax) returns the
%   standard bridge envelope (see CONTRACT.md) carrying:
%       .candidates  --  algorithm-defined result struct (rows/targets/config)
%       .summary     --  short scalar struct with at least a 'count' field
%
%   The bridge defers to ``algOptknock`` (in StrainDesign/algorithms) and
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
    if nargin < 4 || isempty(max_candidates), max_candidates = 200;  end
    if nargin < 5 || isempty(num_del),         num_del         = 5;    end
    if nargin < 6 || isempty(min_growth_fraction), min_growth_fraction = 0.1; end
    if nargin < 7 || isempty(vmax),            vmax            = 1000; end

    try
        model = resolve_model_id(model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end

    bridge_log('mdpRunOptknock', 'target=%s bio=%s maxCand=%d k=%d minMu=%.2f', ...
               char(target), char(biomass), double(max_candidates), ...
               double(num_del), double(min_growth_fraction));

    try
        if exist('algOptknock', 'file') ~= 2 && exist('algOptknock', 'file') ~= 5
            error('algOptknock:notFound', 'algOptknock helper not on path');
        end
        candidates = algOptknock(model, char(biomass), char(target), [], ...
                                 double(num_del), double(max_candidates), ...
                                 double(min_growth_fraction), 'GE', false, '');
    catch err
        result = make_err('err_gurobi_license', err.message);
        return;
    end
    if isempty(candidates) || ~isstruct(candidates)
        candidates = struct();
        count = 0;
    else
        count = 0;
        if isfield(candidates, 'rows') && ~isempty(candidates.rows)
            count = numel(candidates.rows);
        elseif isfield(candidates, 'targets') && ~isempty(candidates.targets)
            count = numel(candidates.targets);
        end
    end
    payload = struct('candidates', candidates, ...
                     'summary',    struct('count', count, ...
                                         'max_candidates', double(max_candidates), ...
                                         'num_del', double(num_del), ...
                                         'min_growth_fraction', double(min_growth_fraction), ...
                                         'vmax', double(vmax)));
    result = make_ok(payload);
end
