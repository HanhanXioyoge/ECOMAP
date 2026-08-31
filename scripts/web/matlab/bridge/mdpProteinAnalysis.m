function result = mdpProteinAnalysis(ec_model_id, group_by)
%MDPPROTEINANALYSIS Aggregate protein-usage statistics by subsystem.
%   result = mdpProteinAnalysis(ec_model_id, group_by) returns the
%   standard bridge envelope (see CONTRACT.md) carrying:
%       .usage_table  --  array of structs {subsystem, share} per group
%       .pie_data     --  array of structs {label, value}   suitable for a pie
%
%   `group_by` is one of 'subsystem' (default), 'compartment', or 'gene';
%   unknown values fall back to 'subsystem'. The share per group is the
%   sum of kcat weights (proxy for enzyme demand) divided by the total.
    here = fileparts(mfilename('fullpath'));
    addpath_once(here);                                                % bridge/ for _bridge_util
    addpath_once(fullfile(here, '..', '..', 'Analysis'));              % scripts/Analysis

    if nargin < 1 || isempty(ec_model_id)
        result = make_err('err_param_invalid', 'ec_model_id required');
        return;
    end

    try
        ecModel = resolve_model_id(ec_model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end

    if nargin < 2 || isempty(group_by)
        group_by = 'subsystem';
    end
    if ~ischar(group_by) && ~(isstring(group_by) && isscalar(group_by))
        group_by = 'subsystem';
    else
        group_by = char(lower(group_by));
    end
    if ~ismember(group_by, {'subsystem', 'compartment', 'gene'})
        group_by = 'subsystem';
    end

    bridge_log('mdpProteinAnalysis', 'group_by=%s', group_by);

    if ~isfield(ecModel, 'rxns') || isempty(ecModel.rxns) || ...
       ~isfield(ecModel, 'enzymeConstraints') || ...
       ~isfield(ecModel.enzymeConstraints, 'kcat') || ...
       isempty(ecModel.enzymeConstraints.kcat)
        payload = struct('usage_table', struct('subsystem', {}, 'share', {}), ...
                         'pie_data',    struct('label', {}, 'value', {}));
        result = make_ok(payload);
        return;
    end

    nRxns = numel(ecModel.rxns);
    weight = ecModel.enzymeConstraints.kcat;
    if numel(weight) < nRxns
        weight(end+1:nRxns) = 0;
    end
    weight = weight(1:nRxns);

    switch group_by
        case 'compartment'
            groups = repmat({'cytosol'}, nRxns, 1);
            if isfield(ecModel, 'compartments') && ~isempty(ecModel.compartments)
                for k = 1:min(nRxns, numel(ecModel.compartments))
                    c = ecModel.compartments{k};
                    if iscell(c) && ~isempty(c), c = c{1}; end
                    if ischar(c), groups{k} = char(c); end
                end
            end
        case 'gene'
            groups = repmat({'unassigned'}, nRxns, 1);
            if isfield(ecModel, 'genes') && ~isempty(ecModel.genes)
                for k = 1:nRxns
                    if iscell(ecModel.rxns), rx = char(ecModel.rxns{k});
                    else, rx = sprintf('rxn_%d', k); end
                    rx = regexprep(rx, '(_REV)?(_EXP_\d+)?$', '');
                    gpos = findRxnIDs(ecModel, rx);
                    if isempty(gpos), gpos = k; end
                    if gpos <= numel(ecModel.genes)
                        if iscell(ecModel.genes)
                            groups{k} = char(ecModel.genes{gpos});
                        else
                            groups{k} = sprintf('gene_%d', gpos);
                        end
                    end
                end
            end
        otherwise  % 'subsystem'
            groups = repmat({'unassigned'}, nRxns, 1);
            if isfield(ecModel, 'subSystems') && ~isempty(ecModel.subSystems)
                for k = 1:min(nRxns, numel(ecModel.subSystems))
                    s = ecModel.subSystems{k};
                    if iscell(s) && ~isempty(s), s = s{1}; end
                    if ischar(s), groups{k} = char(s); end
                end
            end
    end
    total = sum(weight);
    if total < eps, total = 1; end
    keys = unique(groups);
    table_out = struct('subsystem', {}, 'share', {});
    pie_out   = struct('label', {}, 'value', {});
    for k = 1:numel(keys)
        mask = strcmp(groups, keys{k});
        share = sum(weight(mask)) / total;
        table_out(k).subsystem = char(keys{k});
        table_out(k).share     = share;
        pie_out(k).label = char(keys{k});
        pie_out(k).value = share;
    end
    result = make_ok(struct('usage_table', table_out, 'pie_data', pie_out));
end
