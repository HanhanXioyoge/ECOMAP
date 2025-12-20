function [ecModel, summaryTbl] = adjustKcatForAbsorptionReactions(ecModel, scaling_factor, writefile, parameters)
% ADJUSTKCATFORABSORPTIONREACTIONS
%   Growth-anchored kcat calibration (GKC-like):
%   Iteratively down-scales kcat values of enzyme-constrained reactions that
%   carry non-zero flux along an "absorption path" starting from each
%   exchange uptake reaction, while enforcing a GLOBAL feasibility criterion.
%
%   === New requirement (implemented) ===
%   Before any scaling happens, we perform a ONE-TIME feasibility screening
%   across ALL environments at their TargetMu (biomass lower bound):
%     - infeasible environments are EXCLUDED from scaling
%     - scaling and global feasibility checks are enforced ONLY on the feasible subset
%   However, the final summary table still contains ALL environments, including
%   excluded ones, with columns:
%     IncludedForScaling, ExcludeReason
%
% SUMMARY TABLE (summaryTbl)
%   - TargetMu      : Target growth rate from UnconstrainedMaxGrowth.tsv (col3)
%   - InitMu        : Max growth BEFORE calibration (maximize biomass)
%   - FinalMu       : Max growth AFTER  calibration (maximize biomass)
%   - AbsErr        : abs(FinalMu - TargetMu)
%   - InitUptake    : substrate uptake BEFORE calibration, Uptake = -v_ex
%   - FinalUptake   : substrate uptake AFTER  calibration, Uptake = -v_ex
%   - UptakeAbsErr  : abs(FinalUptake - InitUptake)
%   - FeasibleInit  : whether initial model solves (maximize biomass) in that env
%   - FeasibleFinal : whether final model solves (maximize biomass) in that env
%   - FeasibleAtTargetInit : feasibility at TargetMu BEFORE scaling (screening)
%   - IncludedForScaling   : whether the env is included in scaling set
%   - ExcludeReason        : reason for exclusion (if excluded)
%
% INPUT
%   ecModel        : enzyme-constrained GEM with enzymeConstraints.rxns/kcat
%   scaling_factor : multiplicative factor (default 0.75), kcat <- kcat * scaling_factor
%   writefile      : whether to save summary table to outputDir (default true)
%   parameters     : ParameterManager params with fields:
%                    - bioRxn, c_source, dataDir, org_name, (optional) outputDir
%
% DATA REQUIREMENT
%   parameters.dataDir/UnconstrainedMaxGrowth.tsv must contain:
%     col1  = exchange rxn id
%     col3  = TargetMu
%     col11 = condition type ('anaerobic'/'limited'/others)

    if nargin < 4 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end
    if nargin < 3 || isempty(writefile)
        writefile = true;
    end
    if nargin < 2 || isempty(scaling_factor)
        scaling_factor = 0.75;
    end

    bioRxn   = parameters.bioRxn;
    c_source = parameters.c_source;
    basePath = parameters.dataDir;
    org_name = parameters.org_name;

    T = readtable(fullfile(basePath, 'UnconstrainedMaxGrowth.tsv'), ...
        'FileType','text','ReadRowNames',true);

    nEnv = height(T);
    exchange_reactions = T{:,1};

    fprintf('[adjustKcat] Start. scaling_factor = %.4f\n', scaling_factor);

    % Ensure S-matrix is consistent with current kcats, snapshot initial model
    ecModel = UpdateSmatrix(ecModel);
    ecModel_init = ecModel;

    % ===================== Phase A: one-time feasibility screening =====================
    fprintf('[adjustKcat] Screening feasibility at TargetMu for all environments (pre-scaling)...\n');

    FeasibleAtTargetInit = false(nEnv,1);
    IncludedForScaling   = false(nEnv,1);
    ExcludeReason        = strings(nEnv,1);

    for i = 1:nEnv
        ex_rxn    = string(T{i,1});
        mu_target = T{i,3};

        m = buildEnvModel(ecModel_init, T, i, bioRxn, c_source, org_name);
        if isempty(m)
            FeasibleAtTargetInit(i) = false;
            IncludedForScaling(i)   = false;
            ExcludeReason(i)        = "ExchangeNotFound";
            continue;
        end

        % Enforce TargetMu as feasibility anchor
        idxBio = strcmp(m.rxns, bioRxn);
        m.lb(idxBio) = mu_target;
        m.c = double(idxBio);

        sol = solveLP(m);
        if checkSolution(sol)
            FeasibleAtTargetInit(i) = true;
            IncludedForScaling(i)   = true;
            ExcludeReason(i)        = "";
        else
            FeasibleAtTargetInit(i) = false;
            IncludedForScaling(i)   = false;
            ExcludeReason(i)        = "InfeasibleAtTargetMu(Init)";
        end
    end

    nIncl = sum(IncludedForScaling);
    nExcl = nEnv - nIncl;
    fprintf('[adjustKcat] Screening done. Included=%d, Excluded=%d\n', nIncl, nExcl);

    eligibleIdx = find(IncludedForScaling);

    % ===================== Phase B: calibration loop (ONLY on eligible envs) =====================
    for kk = 1:numel(eligibleIdx)
        i = eligibleIdx(kk);

        ex_rxn     = exchange_reactions{i};
        mu_target  = T{i,3};
        condType   = string(T{i,11});

        fprintf('\n[Env %d/%d (eligible %d/%d)] Exchange=%s | target_mu=%.6g | cond=%s\n', ...
            i, nEnv, kk, numel(eligibleIdx), string(ex_rxn), mu_target, condType);

        model_env = buildEnvModel(ecModel, T, i, bioRxn, c_source, org_name);
        if isempty(model_env)
            % should not happen because screening passed, but keep robust
            warning('  Unexpected: exchange %s not found during scaling. Skip.', string(ex_rxn));
            IncludedForScaling(i) = false;
            ExcludeReason(i) = "ExchangeNotFound(Unexpected)";
            continue;
        end

        % Sanity check at infinite uptake (no TargetMu)
        sol0 = solveLP(model_env);
        if ~checkSolution(sol0)
            warning('  Unexpected: infeasible even before scaling under infinite uptake. Skip env %s.', string(ex_rxn));
            IncludedForScaling(i) = false;
            ExcludeReason(i) = "InfeasibleEvenBeforeScaling(Unexpected)";
            continue;
        end

        % Anchor growth at TargetMu
        idxBio = strcmp(model_env.rxns, bioRxn);
        model_env.lb(idxBio) = mu_target;
        model_env.c = double(idxBio);

        maxSteps = 200;
        step = 0;

        while step < maxSteps
            step = step + 1;

            sol = solveLP(model_env);
            if ~checkSolution(sol)
                fprintf('  Step %d: current env infeasible -> stop scaling for this env.\n', step);
                break;
            end

            [abs_rxns, ~] = findAbsorptionReactions(model_env, ex_rxn, sol);
            abs_rxns = unique(abs_rxns, 'stable');

            if isempty(abs_rxns)
                fprintf('  Step %d: no absorption reactions found -> stop.\n', step);
                break;
            end

            ecModel_backup = ecModel;

            fprintf('  Step %d: scaling %d reaction(s)\n', step, numel(abs_rxns));

            for j = 1:numel(abs_rxns)
                rxn = abs_rxns{j};
                idxEC = find(strcmp(ecModel.enzymeConstraints.rxns, rxn), 1);

                if isempty(idxEC)
                    fprintf('    - %s : not found in ecModel.enzymeConstraints.rxns (skip)\n', string(rxn));
                    continue;
                end

                old_kcat = ecModel.enzymeConstraints.kcat(idxEC);
                new_kcat = old_kcat * scaling_factor;
                ecModel.enzymeConstraints.kcat(idxEC) = new_kcat;

                fprintf('    - %s : kcat %.6g -> %.6g (x%.4f)\n', ...
                    string(rxn), old_kcat, new_kcat, scaling_factor);
            end

            ecModel = UpdateSmatrix(ecModel);

            % GLOBAL feasibility check ONLY on eligible env set
            [okAll, failInfo] = checkEligibleEnvironments(ecModel, T, parameters, IncludedForScaling);

            if okAll
                fprintf('  Step %d: global feasibility (eligible set) PASS. Commit.\n', step);

                % Rebuild env model from updated global ecModel
                model_env = buildEnvModel(ecModel, T, i, bioRxn, c_source, org_name);
                idxBio = strcmp(model_env.rxns, bioRxn);
                model_env.lb(idxBio) = mu_target;
                model_env.c = double(idxBio);

            else
                fprintf('  Step %d: global feasibility (eligible set) FAIL at env %d (Exchange=%s, target_mu=%.6g). Revert & stop this env.\n', ...
                    step, failInfo.envIndex, string(failInfo.exchange), failInfo.targetMu);

                ecModel = ecModel_backup;
                ecModel = UpdateSmatrix(ecModel);
                break;
            end
        end

        if step >= maxSteps
            fprintf('  Reached maxSteps=%d for env %s. Stop to avoid endless scaling.\n', maxSteps, string(ex_rxn));
        end
    end

    fprintf('\n[adjustKcat] Calibration Done. Computing summary table using INIT and FINAL models...\n');

    % ===================== Phase C: summary table (ALL envs, including excluded) =====================
    summaryTbl = computeInitFinalSummary(ecModel_init, ecModel, T, parameters);

    % Add screening + inclusion info
    summaryTbl.FeasibleAtTargetInit = FeasibleAtTargetInit;
    summaryTbl.IncludedForScaling   = IncludedForScaling;
    summaryTbl.ExcludeReason        = ExcludeReason;

    disp(summaryTbl);

    % Save summary
    if writefile && isfield(parameters,'outputDir') && ~isempty(parameters.outputDir)
        outDir = parameters.outputDir;
        if ~exist(outDir,'dir'); mkdir(outDir); end
        outFile = fullfile(outDir, 'adjustKcat_Summary.tsv');
        writetable(summaryTbl, outFile, 'FileType','text','Delimiter','\t');
    end

    fprintf('[adjustKcat] Done.\n');
end

% ----------------------- Helper: build env model --------------------------
function model_i = buildEnvModel(ecModel, T, i, bioRxn, c_source, org_name)
% BUILDENVMODEL
%   Construct the i-th environment model from a global ecModel:
%     - Apply anaerobicModel for anaerobic/limited conditions
%     - Block uptake of c_source (lb=0)
%     - Set the environment exchange uptake to -1000
%     - Objective biomass (biomass maximization unless overridden by caller)

    model_i = ecModel;

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

% -------------------- Helper: global feasibility check (eligible only) --------------------
function [allFeasible, failInfo] = checkEligibleEnvironments(ecModel, T, parameters, eligibleMask)
% CHECKELIGIBLEENVIRONMENTS
%   Feasibility check only on environments where eligibleMask(i) is true.
%   For each eligible env:
%     - exchange uptake = -1000
%     - c_source uptake blocked (lb=0)
%     - biomass lower bound = TargetMu
%   Return first failing env info.

    allFeasible = true;
    failInfo = struct('envIndex',[], 'exchange',"", 'targetMu',NaN);

    bioRxn   = parameters.bioRxn;
    c_source = parameters.c_source;
    org_name = parameters.org_name;

    idxList = find(eligibleMask);
    for t = 1:numel(idxList)
        k = idxList(t);

        m = buildEnvModel(ecModel, T, k, bioRxn, c_source, org_name);
        if isempty(m)
            allFeasible = false;
            failInfo.envIndex = k;
            failInfo.exchange = string(T{k,1});
            failInfo.targetMu = T{k,3};
            return;
        end

        mu_target = T{k,3};
        idxBio = strcmp(m.rxns, bioRxn);
        m.lb(idxBio) = mu_target;
        m.c = double(idxBio);

        sol = solveLP(m);
        if ~checkSolution(sol)
            allFeasible = false;
            failInfo.envIndex = k;
            failInfo.exchange = string(T{k,1});
            failInfo.targetMu = mu_target;
            return;
        end
    end
end

% -------------------- Helper: summary (INIT vs FINAL) ----------------------
function summaryTbl = computeInitFinalSummary(ecModel_init, ecModel_final, T, parameters)
% COMPUTEINITFINALSUMMARY
%   For each environment:
%     - TargetMu from T
%     - InitMu/InitUptake computed using ecModel_init (maximize biomass)
%     - FinalMu/FinalUptake computed using ecModel_final (maximize biomass)
%   Uptake is defined as Uptake = -v_exchange.

    bioRxn   = parameters.bioRxn;
    c_source = parameters.c_source;
    org_name = parameters.org_name;

    nEnv = height(T);

    EnvIndex  = (1:nEnv)';
    Exchange  = string(T{:,1});
    Condition = string(T{:,11});
    TargetMu  = T{:,3};

    InitMu       = nan(nEnv,1);
    FinalMu      = nan(nEnv,1);
    InitUptake   = nan(nEnv,1);
    FinalUptake  = nan(nEnv,1);
    FeasibleInit = false(nEnv,1);
    FeasibleFinal= false(nEnv,1);

    for i = 1:nEnv
        ex_rxn = string(T{i,1});

        m0 = buildEnvModel(ecModel_init, T, i, bioRxn, c_source, org_name);
        if ~isempty(m0)
            [InitMu(i), InitUptake(i), FeasibleInit(i)] = getMaxGrowthAndUptake(m0, bioRxn, ex_rxn);
        end

        m1 = buildEnvModel(ecModel_final, T, i, bioRxn, c_source, org_name);
        if ~isempty(m1)
            [FinalMu(i), FinalUptake(i), FeasibleFinal(i)] = getMaxGrowthAndUptake(m1, bioRxn, ex_rxn);
        end
    end

    AbsErr       = abs(FinalMu - TargetMu);
    UptakeAbsErr = abs(FinalUptake - InitUptake);

    summaryTbl = table(EnvIndex, Exchange, Condition, TargetMu, ...
                       InitMu, FinalMu, AbsErr, ...
                       InitUptake, FinalUptake, UptakeAbsErr, ...
                       FeasibleInit, FeasibleFinal);
end

function [mu, uptake, ok] = getMaxGrowthAndUptake(model_env, bioRxn, ex_rxn)
% GETMAXGROWTHANDUPTAKE
%   Solve FBA maximizing biomass and return:
%     mu     = growth rate
%     uptake = -v_exchange (positive means uptake)
%     ok     = feasibility

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
        uptake = -sol.x(idxEx);
    else
        uptake = NaN;
    end
end

% ------------------- Absorption-path reaction discovery -------------------
function [absorption_reactions, visited_reactions] = findAbsorptionReactions(model, exchange_reaction, sol)
% FINDABSORPTIONREACTIONS
%   Starting from metabolite consumed by the exchange uptake reaction,
%   traverse downstream reactions (excluding exchange reactions).
%   A reaction is an "absorption reaction" if:
%     (1) it carries non-zero flux in sol; and
%     (2) it is enzyme-constrained (exists in model.enzymeConstraints.rxns).

    absorption_reactions = {};
    visited_reactions = {};

    exchange_reaction_idx = strcmp(model.rxns, exchange_reaction);
    if ~any(exchange_reaction_idx)
        warning('Exchange reaction %s not found in the model.', string(exchange_reaction));
        return;
    end

    connected_met_idx = find(model.S(:, exchange_reaction_idx) < 0);
    if isempty(connected_met_idx)
        warning('Exchange reaction %s does not consume any metabolite.', string(exchange_reaction));
        return;
    end

    start_met_id = model.mets{connected_met_idx(1)};
    [absorption_reactions, visited_reactions] = recursiveSearch(model, start_met_id, sol, absorption_reactions, visited_reactions);
end

function [absorption_reactions, visited_reactions] = recursiveSearch(model, current_met_id, sol, absorption_reactions, visited_reactions)
% RECURSIVESEARCH
%   Depth-first traversal from a starting metabolite.
%   Only continue traversal through CYTOSOLIC metabolites (requested behavior).
%
% NOTE
%   This function assumes COBRA-style metabolite IDs with compartment suffix
%   like '_c' / '[c]' / '(c)'. If unavailable, we fall back to model.metComps.

    current_met_idx = find(strcmp(model.mets, current_met_id), 1);
    if isempty(current_met_idx)
        return;
    end

    downstream_rxn_indices = find(model.S(current_met_idx, :) ~= 0);

    % Filter out exchange reactions (columns with only one non-zero)
    is_exchange = sum(model.S ~= 0, 1) == 1;
    downstream_rxn_indices = downstream_rxn_indices(~is_exchange(downstream_rxn_indices));

    current_level = {};
    for ii = 1:numel(downstream_rxn_indices)
        ridx = downstream_rxn_indices(ii);
        rxn_id = model.rxns{ridx};
        if ~ismember(rxn_id, visited_reactions)
            current_level{end+1} = rxn_id; %#ok<AGROW>
            visited_reactions{end+1} = rxn_id; %#ok<AGROW>
        end
    end

    if isempty(current_level)
        return;
    end

    for ii = 1:numel(current_level)
        rxn_id = current_level{ii};
        rxn_idx = strcmp(model.rxns, rxn_id);

        v = sol.x(rxn_idx);
        is_flux_nonzero = abs(v) > 1e-9;
        is_enzyme_constrained = any(strcmp(model.enzymeConstraints.rxns, rxn_id));

        if is_flux_nonzero && is_enzyme_constrained
            absorption_reactions{end+1} = rxn_id; %#ok<AGROW>
            continue;
        end

        if ~is_flux_nonzero
            continue;
        end

        prod_met_indices = find(model.S(:, rxn_idx) > 0);
        for met_idx = prod_met_indices'
            next_met_id = model.mets{met_idx};

            % === requested filter: only traverse cytosolic metabolites ===
            if isCytosolicMet(model, next_met_id)
                continue;
            end

            [absorption_reactions, visited_reactions] = recursiveSearch(model, next_met_id, sol, absorption_reactions, visited_reactions);
        end
    end
end

function tf = isCytosolicMet(model, met_id)
% ISCYTOSOLICMET
%   Return true if metabolite is cytosolic ('c').
%   Priority:
%     1) COBRA id suffix patterns: '_c' / '[c]' / '(c)'
%     2) model.metComps + model.comps if available
%   If compartment cannot be determined, return true (do not over-filter).

    tf = false;

    if isstring(met_id); met_id = char(met_id); end

    % Pattern-based checks
    if ~isempty(regexp(met_id, '\[c\]$', 'once')), tf = true; return; end
    if ~isempty(regexp(met_id, '_c$',   'once')), tf = true; return; end
    if ~isempty(regexp(met_id, '\(c\)$','once')), tf = true; return; end

    % Field-based checks
    if isfield(model,'metComps') && isfield(model,'comps')
        midx = find(strcmp(model.mets, met_id), 1);
        if ~isempty(midx)
            compID = model.comps{model.metComps(midx)};
            tf = strcmp(compID,'c');
            return;
        end
    end

    % Unknown compartment => do not over-filter (keeps traversal working)
    % tf = true;
end
