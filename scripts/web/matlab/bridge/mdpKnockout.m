function result = mdpKnockout(ec_model_id, gene_list, c_source)
%MDPKNOCKOUT Run single-gene knockout enumeration on an ecModel.
%   result = mdpKnockout(ec_model_id, gene_list, c_source) returns the
%   standard bridge envelope (see CONTRACT.md) carrying:
%       .knockout_table  --  array of structs {gene, growth_ratio} per gene
%       .heatmap         --  same array, reused by the front-end heatmap
%
%   `gene_list` may be 'all' (use every gene in the model), a string 'all',
%   the literal word 'all', or a cell/numeric array of gene IDs. The
%   growth ratio is the absolute biomass flux after deleting each gene,
%   divided by the wild-type biomass flux.
    here = fileparts(mfilename('fullpath'));
    addpath_once(here);                                                % bridge/ for _bridge_util
    addpath_once(fullfile(here, '..', '..', 'Analysis'));              % scripts/Analysis
    addpath_once(fullfile(here, '..', '..', 'StrainDesign'));          % scripts/StrainDesign

    if nargin < 1 || isempty(ec_model_id)
        result = make_err('err_param_invalid', 'ec_model_id required');
        return;
    end
    if nargin < 3 || isempty(c_source)
        result = make_err('err_param_invalid', 'c_source required');
        return;
    end

    try
        ecModel = resolve_model_id(ec_model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end

    if ischar(gene_list) || (isstring(gene_list) && isscalar(gene_list))
        gl = char(gene_list);
        if strcmpi(gl, 'all')
            genes = ecModel.genes;
        else
            genes = {gl};
        end
    elseif iscell(gene_list)
        genes = gene_list;
    elseif isnumeric(gene_list)
        if gene_list == 0, genes = ecModel.genes; else, genes = ecModel.genes(gene_list); end
    else
        genes = ecModel.genes;
    end
    if isempty(genes)
        result = make_err('err_param_invalid', 'no genes to knock out');
        return;
    end
    bridge_log('mdpKnockout', 'genes=%d source=%s', numel(genes), char(c_source));

    try
        wt = setParam(ecModel, 'lb', c_source, -10);
        wtSol = solveLP(wt);
    catch err
        result = make_err('err_gurobi_license', err.message);
        return;
    end
    wtGrowth = abs(wtSol.f);
    if wtGrowth < eps, wtGrowth = 1; end

    n = numel(genes);
    table = struct('gene', {}, 'growth_ratio', {});
    for k = 1:n
        if iscell(genes), g = char(genes{k}); else, g = sprintf('gene_%d', k); end
        try
            m = ecModel;
            m.lb(findRxnIDs(m, c_source)) = -10;
            rxnIdxs = findRxnIDs(m, g);
            m.lb(rxnIdxs) = 0;
            m.ub(rxnIdxs) = 0;
            sol = solveLP(m);
            if ~isempty(sol) && isfield(sol, 'f') && ~isempty(sol.f) && abs(sol.f) >= 0
                ratio = abs(sol.f) / wtGrowth;
                table(k).gene = g;
                table(k).growth_ratio = ratio;
            end
        catch errKO
            bridge_log('mdpKnockout', 'KO %s failed: %s', g, errKO.message);
        end
    end
    result = make_ok(struct('knockout_table', table, 'heatmap', table));
end
