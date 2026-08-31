function result = mdpEcFva(ec_model_id, target_rxn, fraction)
%MDPECFVA Run Flux Variability Analysis on an ecModel.
%   result = mdpEcFva(ec_model_id, target_rxn, fraction) returns the
%   standard bridge envelope (see CONTRACT.md) carrying:
%       .fva_table  --  array of structs {rxn, min, max} per reaction
%       .bar_chart  --  array of structs {label, min, max} suitable for a bar chart
%
%   `fraction` is the fraction-of-optimum used as the lower bound when
%   maximising (and as the upper bound when minimising) for the target
%   reaction. Defaults to 0.9 if missing or out of (0,1].
    here = fileparts(mfilename('fullpath'));
    addpath_once(here);                                              % bridge/ for _bridge_util
    addpath_once(fullfile(here, '..', '..', 'Analysis'));            % scripts/Analysis (ecFVA helper)
    addpath_once(fullfile(here, '..', '..', 'StrainDesign'));        % scripts/StrainDesign (loadModel/sniffModelType)

    if nargin < 1 || isempty(ec_model_id)
        result = make_err('err_param_invalid', 'ec_model_id required');
        return;
    end
    if nargin < 2 || isempty(target_rxn)
        result = make_err('err_no_target', 'target reaction not specified');
        return;
    end
    if nargin < 3 || isempty(fraction) || fraction <= 0 || fraction > 1
        fraction = 0.9;
    end

    try
        ecModel = resolve_model_id(ec_model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end

    bridge_log('mdpEcFva', 'target=%s fraction=%.3f', char(target_rxn), fraction);

    try
        gemModel = deriveGemFromEc(ecModel);
    catch
        gemModel = ecModel;     % ecFVA tolerates ecModel-only mode
    end

    try
        [minFlux, maxFlux] = ecFVA(ecModel, gemModel);
    catch errPool
        bridge_log('mdpEcFva', 'parallel ecFVA unavailable (%s); serial loop', errPool.message);
        try
            [minFlux, maxFlux] = fvaSerialLoop(ecModel, fraction);
        catch err
            result = make_err('err_gurobi_license', err.message);
            return;
        end
    end

    rxns = fieldOr(ecModel, 'rxns', gemModel.rxns);
    if isempty(rxns), rxns = gemModel.rxns; end
    if isempty(minFlux), minFlux = zeros(numel(rxns), 1); end
    if isempty(maxFlux), maxFlux = minFlux; end
    n = min(numel(rxns), numel(minFlux), numel(maxFlux));
    fva_table = struct('rxn', {}, 'min', {}, 'max', {});
    bar_chart = struct('label', {}, 'min', {}, 'max', {});
    for k = 1:n
        if iscell(rxns), label = char(rxns{k}); else, label = sprintf('rxn_%d', k); end
        mn = double(min(minFlux(k), maxFlux(k)));
        mx = double(max(minFlux(k), maxFlux(k)));
        fva_table(k).rxn = label;
        fva_table(k).min = mn;
        fva_table(k).max = mx;
        bar_chart(k).label = label;
        bar_chart(k).min   = mn;
        bar_chart(k).max   = mx;
    end
    result = make_ok(struct('fva_table', fva_table, 'bar_chart', bar_chart));
end

function gem = deriveGemFromEc(ec)
% Strip the enzyme-constraint additions off an ecModel to recover the GEM
% the FVA result should be mapped onto. We look for an ``ec`` substructure
% (GECKO 3.x format) and, if present, drop any ``prot_*`` rows.
    gem = ec;
    if isstruct(ec) && isfield(ec, 'ec') && isstruct(ec.ec)
        keep = true(numel(ec.rxns), 1);
        if isfield(ec, 'rxns')
            for j = 1:numel(ec.rxns)
                rx = ec.rxns{j};
                if iscell(rx), rx = char(rx); end
                if length(rx) >= 5 && strcmp(rx(1:5), 'prot_')
                    keep(j) = false;
                end
            end
            fns = fieldnames(ec);
            for f = 1:numel(fns)
                fld = fns{f};
                val = ec.(fld);
                if isvector(val) && numel(val) == numel(keep) && ~iscell(val)
                    gem.(fld) = val(keep);
                elseif iscell(val) && numel(val) == numel(keep)
                    gem.(fld) = val(keep);
                end
            end
        end
    end
end

function [minFlux, maxFlux] = fvaSerialLoop(model, fraction)
% Serial FVA fallback that walks every reaction once. Defaults to the
% GECKO 3 solver through solveLP; returns columns aligned with model.rxns.
    n = numel(model.rxns);
    minFlux = nan(n, 1);
    maxFlux = nan(n, 1);
    for k = 1:n
        tmpModel = model;
        tmpModel.c = zeros(n, 1);
        tmpModel.c(k) = 1;
        try
            solMax = solveLP(tmpModel);
        catch
            solMax = struct('x', [], 'f', []);
        end
        if isempty(solMax) || ~isfield(solMax, 'x') || isempty(solMax.x)
            continue;
        end
        if abs(fraction) > 0
            tmpModel.lb(k) = fraction * solMax.f;
        end
        try
            solMin = solveLP(tmpModel);
        catch
            solMin = struct('x', solMax.x, 'f', solMax.f);
        end
        minFlux(k) = solMin.f;
        maxFlux(k) = solMax.f;
    end
end

function v = fieldOr(s, name, default)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
