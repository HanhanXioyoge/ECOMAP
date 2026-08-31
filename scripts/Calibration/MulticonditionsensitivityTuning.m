function [model, summaryTbl, changesTbl] = MulticonditionsensitivityTuning(model, foldChange, parameters)
% MULTICONDITIONSENSITIVITYTUNING
%   Run sensitivity-driven kcat tuning sequentially across multiple environments
%   defined in UnconstrainedMaxGrowth.tsv, and return:
%     (1) summaryTbl: compares experimental target growth (TargetMu) with:
%         - InitMu      : predicted maximal growth BEFORE tuning
%         - FinalMu     : predicted maximal growth AFTER all tuning ends
%         - AbsErr      : abs(FinalMu - TargetMu)
%         - InitUptake  : predicted substrate uptake rate BEFORE tuning (Uptake=-v_ex)
%         - FinalUptake : predicted substrate uptake rate AFTER tuning  (Uptake=-v_ex)
%         - UptakeAbsErr: abs(FinalUptake - InitUptake)
%     (2) changesTbl: aggregated change log of all kcat modifications, labeled
%         by EnvIndex / Exchange / Condition / Iteration (if available).
%
%   Important design choice:
%     Because later tuning steps may change kcats and thus affect earlier environments,
%     we DO NOT report per-environment MuAfter right after its own tuning step.
%     Instead, we compute InitMu/InitUptake using the initial model and
%     FinalMu/FinalUptake using the final model, both evaluated across ALL environments.
%
% INPUT
%   model      : ecModel with fields enzymeConstraints.kcat and enzymeConstraints.rxns
%   foldChange : multiplicative factor used in sensitivityTuning (default 10)
%   parameters : ParameterManager params struct with fields:
%                - bioRxn, c_source, dataDir, org_name, (optional) outputDir
%
% OUTPUT
%   model      : final tuned global model (kcat updated)
%   summaryTbl : EnvIndex, Exchange, Condition, TargetMu, InitMu, FinalMu,
%                AbsErr, InitUptake, FinalUptake, UptakeAbsErr,
%                FeasibleInit, FeasibleFinal, NumKcatChangesThisEnv
%   changesTbl : aggregated kcat change log across all environments

    if nargin < 3 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end

    if nargin < 2 || isempty(foldChange)
        foldChange = 10;
    end

    bioRxn   = parameters.bioRxn;
    c_source = parameters.c_source;
    basePath = parameters.calibrationDir;
    org_name = parameters.org_name;

    T = readtable(fullfile(basePath, 'data','UnconstrainedMaxGrowth.tsv'), ...
        'FileType', 'text', 'ReadRowNames', true);

    nEnv = height(T);

    % Ensure objective is biomass on the global model
    idxBio_global = strcmp(model.rxns, bioRxn);
    model.c = double(idxBio_global);

    % Make sure S-matrix matches current kcats
    model = UpdateSmatrix(model);

    % =================== Phase 0: compute InitMu & InitUptake using the INITIAL model ===================
    model_init = model;  % snapshot BEFORE any tuning

    InitMu       = nan(nEnv,1);
    InitUptake   = nan(nEnv,1);
    FeasibleInit = false(nEnv,1);

    for i = 1:nEnv
        ex_rxn = string(T{i,1});

        model_env_init = buildEnvModel(model_init, T, i, bioRxn, c_source, org_name);
        if isempty(model_env_init)
            InitMu(i) = NaN;
            InitUptake(i) = NaN;
            FeasibleInit(i) = false;
        else
            [InitMu(i), InitUptake(i), FeasibleInit(i)] = getMaxGrowthAndUptake(model_env_init, bioRxn, ex_rxn);
        end
    end

    % =================== Phase 1: sequential tuning across environments ===================
    changesTbl = table();
    numChangesThisEnv = zeros(nEnv,1);

    for i = 1:nEnv
        ex_rxn    = string(T{i,1});
        mu_target = T{i,3};
        condType  = string(T{i,11});

        % Build environment model from CURRENT global model
        model_env = buildEnvModel(model, T, i, bioRxn, c_source, org_name);
        if isempty(model_env)
            numChangesThisEnv(i) = 0;
            continue;
        end

        % Backup kcats before tuning this environment (for change logging)
        if isfield(model_env,'enzymeConstraints') && isfield(model_env.enzymeConstraints,'kcat')
            kcat_before = model_env.enzymeConstraints.kcat;
        else
            kcat_before = [];
        end

        % Run tuning to reach target growth
        try
            [model_env_tuned, tunedKcats] = sensitivityTuning(model_env, mu_target, foldChange);
        catch
            numChangesThisEnv(i) = 0;
            continue;
        end

        % kcat after tuning this environment
        if isfield(model_env_tuned,'enzymeConstraints') && isfield(model_env_tuned.enzymeConstraints,'kcat')
            kcat_after = model_env_tuned.enzymeConstraints.kcat;
        else
            kcat_after = [];
        end

        % Count changes for this env tuning step (local diff)
        if ~isempty(kcat_before) && ~isempty(kcat_after) && numel(kcat_before)==numel(kcat_after)
            numChangesThisEnv(i) = sum(abs(kcat_after(:) - kcat_before(:)) > 0);
        else
            numChangesThisEnv(i) = 0;
        end

        % Build change log table for this environment and append
        if ~isempty(kcat_before) && ~isempty(kcat_after) && isfield(model_env_tuned,'enzymeConstraints') ...
                && isfield(model_env_tuned.enzymeConstraints,'rxns')
            envChangeTbl = buildChangeTable( ...
                tunedKcats, model_env_tuned.enzymeConstraints.rxns, ...
                kcat_before, kcat_after, ...
                i, ex_rxn, condType);
            if ~isempty(envChangeTbl)
                changesTbl = [changesTbl; envChangeTbl]; %#ok<AGROW>
            end
        end

        % Sync tuned kcats back to GLOBAL model, rebuild S-matrix
        if isfield(model_env_tuned,'enzymeConstraints') && isfield(model_env_tuned.enzymeConstraints,'kcat')
            model.enzymeConstraints.kcat = model_env_tuned.enzymeConstraints.kcat;
            model = UpdateSmatrix(model);
        end
    end

    % =================== Phase 2: compute FinalMu & FinalUptake using FINAL model across all envs ===================
    FinalMu       = nan(nEnv,1);
    FinalUptake   = nan(nEnv,1);
    FeasibleFinal = false(nEnv,1);

    for i = 1:nEnv
        ex_rxn = string(T{i,1});

        model_env_final = buildEnvModel(model, T, i, bioRxn, c_source, org_name);
        if isempty(model_env_final)
            FinalMu(i) = NaN;
            FinalUptake(i) = NaN;
            FeasibleFinal(i) = false;
        else
            [FinalMu(i), FinalUptake(i), FeasibleFinal(i)] = getMaxGrowthAndUptake(model_env_final, bioRxn, ex_rxn);
        end
    end

    % =================== Build summaryTbl ===================
    EnvIndex  = (1:nEnv)';
    Exchange  = string(T{:,1});
    Condition = string(T{:,11});
    TargetMu  = T{:,3};

    AbsErr        = abs(FinalMu - TargetMu);
    UptakeAbsErr  = abs(FinalUptake - InitUptake);

    summaryTbl = table(EnvIndex, Exchange, Condition, TargetMu, InitMu, FinalMu, AbsErr, ...
                       InitUptake, FinalUptake, UptakeAbsErr, ...
                       FeasibleInit, FeasibleFinal, numChangesThisEnv, ...
                       'VariableNames', {'EnvIndex','Exchange','Condition','TargetMu', ...
                                         'InitMu','FinalMu','AbsErr', ...
                                         'InitUptake','FinalUptake','UptakeAbsErr', ...
                                         'FeasibleInit','FeasibleFinal','NumKcatChangesThisEnv'});

    % ------------------ Only print the two tables ------------------
    disp(summaryTbl);
    disp(changesTbl);

    % Optional: save results to files (silent)
    if isfield(parameters,'outputDir') && ~isempty(parameters.outputDir)
        outDir = parameters.calibrationDir;
        if ~exist(outDir,'dir'); mkdir(outDir); end
        writetable(summaryTbl, fullfile(outDir,'MultiCondition_Summary.tsv'), 'FileType','text','Delimiter','\t');
        writetable(changesTbl, fullfile(outDir,'MultiCondition_KcatChanges.tsv'), 'FileType','text','Delimiter','\t');
    end
end

% ====================== Helper: build environment model ======================
function model_i = buildEnvModel(model, T, i, bioRxn, c_source, org_name)
% BUILDENVMODEL
%   Construct an environment-specific model from a global model:
%     - Apply anaerobic/limited constraints if needed
%     - Block uptake of c_source (lb=0)
%     - Set the environment exchange uptake lower bound to -1000 (unbounded uptake)
%     - Objective is biomass

    model_i = model;

    ex_rxn   = string(T{i,1});
    condType = string(T{i,11});

    if condType=="anaerobic" || condType=="limited"
        model_i = anaerobicModel(model_i, org_name);
    end

    idxC = strcmp(model_i.rxns, c_source);
    if any(idxC), model_i.lb(idxC) = 0; end

    idxEx = strcmp(model_i.rxns, ex_rxn);
    if ~any(idxEx)
        model_i = [];
        return;
    end
    model_i.lb(idxEx) = -1000;

    model_i.c = double(strcmp(model_i.rxns, bioRxn));
end

% ====================== Helper: compute maximal growth + substrate uptake ======================
function [mu, uptake, ok] = getMaxGrowthAndUptake(model_env, bioRxn, ex_rxn)
% GETMAXGROWTHANDUPTAKE
%   Solve FBA maximizing biomass and return:
%     - mu     : objective value (growth rate)
%     - uptake : substrate uptake rate defined as Uptake = -v_ex (COBRA uptake is negative)
%     - ok     : whether solution is feasible

    model_env.c = double(strcmp(model_env.rxns, bioRxn));
    sol = solveLP(model_env);
    ok = checkSolution(sol);

    if ok && isfield(sol,'f') && ~isempty(sol.f)
        mu = sol.f;
    else
        mu = NaN;
    end

    idxEx = strcmp(model_env.rxns, ex_rxn);
    if ok && any(idxEx) && isfield(sol,'x') && ~isempty(sol.x)
        uptake = -sol.x(idxEx);   % uptake as positive if flux is negative
    else
        uptake = NaN;
    end
end

% ====================== Helper: build per-env change table ======================
function envChangeTbl = buildChangeTable(tunedKcats, ecRxns, kcat_before, kcat_after, envIdx, ex_rxn, condType)
% BUILDCHANGETABLE
%   Produce a table of kcat changes for one tuning call (one environment).
%   Priority:
%     1) If tunedKcats contains reaction IDs (and possibly iteration index),
%        use it to generate a compact log with Iteration column.
%     2) Otherwise fallback to diff(kcat) across all ecRxns (Iteration=NaN).
%
% OUTPUT COLUMNS
%   EnvIndex, Exchange, Condition, Iteration, Rxn, OldKcat, NewKcat, FoldChange

    envChangeTbl = table();

    rxnList = {};
    iterList = [];

    if ~isempty(tunedKcats)
        if istable(tunedKcats)
            rxnList  = getFirstExistingColumnAsCellstr(tunedKcats, {'rxn','rxns','reaction','Reaction','Rxn','RxnID'});
            iterList = getFirstExistingColumnAsNumeric(tunedKcats, {'iter','Iter','iteration','Iteration','step','Step'});
        elseif isstruct(tunedKcats)
            try
                TT = struct2table(tunedKcats);
                rxnList  = getFirstExistingColumnAsCellstr(TT, {'rxn','rxns','reaction','Reaction','Rxn','RxnID'});
                iterList = getFirstExistingColumnAsNumeric(TT, {'iter','Iter','iteration','Iteration','step','Step'});
            catch
                rxnList = {};
                iterList = [];
            end
        end
    end

    if ~isempty(rxnList)
        n = numel(rxnList);
        if isempty(iterList), iterList = nan(n,1); end
        if isrow(iterList), iterList = iterList(:); end

        oldK = nan(n,1); newK = nan(n,1);
        for r = 1:n
            idx = find(strcmp(ecRxns, rxnList{r}), 1);
            if ~isempty(idx)
                oldK(r) = kcat_before(idx);
                newK(r) = kcat_after(idx);
            end
        end
        fold = newK ./ oldK;

        envChangeTbl = table( ...
            repmat(envIdx,n,1), repmat(string(ex_rxn),n,1), repmat(string(condType),n,1), ...
            iterList, string(rxnList(:)), oldK, newK, fold, ...
            'VariableNames', {'EnvIndex','Exchange','Condition','Iteration','Rxn','OldKcat','NewKcat','FoldChange'} ...
        );
        return;
    end

    if isempty(kcat_before) || isempty(kcat_after) || numel(kcat_before)~=numel(kcat_after)
        envChangeTbl = table();
        return;
    end

    idxChg = find(abs(kcat_after(:) - kcat_before(:)) > 0);
    if isempty(idxChg)
        envChangeTbl = table();
        return;
    end

    n = numel(idxChg);
    rxnList = string(ecRxns(idxChg));
    oldK = kcat_before(idxChg);
    newK = kcat_after(idxChg);
    fold = newK ./ oldK;

    envChangeTbl = table( ...
        repmat(envIdx,n,1), repmat(string(ex_rxn),n,1), repmat(string(condType),n,1), ...
        nan(n,1), rxnList, oldK, newK, fold, ...
        'VariableNames', {'EnvIndex','Exchange','Condition','Iteration','Rxn','OldKcat','NewKcat','FoldChange'} ...
    );
end

% ====================== Helper: robust column extraction ======================
function vals = getFirstExistingColumnAsCellstr(T, names)
    vals = {};
    for k = 1:numel(names)
        if any(strcmp(T.Properties.VariableNames, names{k}))
            v = T.(names{k});
            if isstring(v), v = cellstr(v); end
            if iscell(v), vals = v; return; end
            if ischar(v), vals = cellstr(v); return; end
        end
    end
end

function vals = getFirstExistingColumnAsNumeric(T, names)
    vals = [];
    for k = 1:numel(names)
        if any(strcmp(T.Properties.VariableNames, names{k}))
            v = T.(names{k});
            if isnumeric(v), vals = v; return; end
        end
    end
end
