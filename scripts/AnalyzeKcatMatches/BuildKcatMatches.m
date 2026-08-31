function MATCH = BuildKcatMatches(DeepLearningModel, complex_name, aggMethod, filePath, parameters)
% BuildKcatMatches
% ----------------
% Purpose
%   Load experimental kcat DB tables (BRENDA + SABIO), load prediction tables
%   (DLKcat / UniKP / CatPred), perform hierarchical matching with a
%   metabolite-name fallback, aggregate duplicate pairs per
%   (ProteinID, metaboliteKey), and return ONLY the matched tables.
%
% Matching logic
%   (1) ProteinID == uniprot
%   (2) sequence identical (only enforced if predicted sequence is non-empty)
%   (3) metabolite identity by ANY of:
%       - exact InChIKey
%       - exact MetaNetXID
%       - exact Substrate_norm (case/spacing normalized name)
%
% Group key (internal)
%   metaboliteKey = IK:<InChIKey> | MNX:<MetaNetXID> | NM:<Substrate_norm>
%   Priority: InChIKey > MetaNetXID > Substrate_norm
%   (Note: this key is used internally for grouping; it is NOT returned.)
%
% Aggregation within each (ProteinID, metaboliteKey) group
%   By default: arithmetic MEAN in log10-space for predicted_kcat_log10 and exp_kcat_log10.
%   Switch to 'median' via aggMethod = 'median'.
%
% Complex flag
%   Extra input 'complex_name' lists reaction names considered "complex".
%   If the aggregated ReactionName equals any name in complex_name (case- /
%   whitespace-insensitive), the output column 'isComplex' is true.
%
% Output
%   MATCH : struct with one field per model tag (DLKcat/UniKP/CatPred).
%           Each field is a table with columns:
%             ReactionName, ProteinID, ec, MetaNetXID, InChIKey, Substrate,
%             n_group,
%             predicted_kcat_log10, exp_kcat_log10,
%             predicted_kcat, exp_kcat,
%             isComplex
%           (No metaboliteKey is returned; Substrate is the human-readable name.)
%
% Inputs
%   DeepLearningModel : string|char|cellstr
%       One or more of {'DLKcat','UniKP','CatPred'}; default = all three.
%   filePath : directory containing <Model>.csv prediction files
%   parameters : struct; if empty, ParameterManager.getParams() is used.
%       - parameters.reconstructionDir (used only if filePath is empty)
%   aggMethod : 'mean' (default) or 'median'  -- aggregation in log10-space
%   complex_name : string|char|cellstr
%       List of reaction names to be flagged as complex. Optional.
%
% Notes
%   - Requires helper functions: findECOMAProot, ParameterManager (your project).
%   - This function ONLY returns matches; no stats/plots/saving.

    % ------------------------- Defaults & args -------------------------
    if nargin < 2, complex_name = []; end
    if nargin < 3 || isempty(aggMethod), aggMethod = 'mean'; end
    if ~any(strcmpi(aggMethod, {'mean','median'}))
        error('BuildKcatMatches:BadArg','aggMethod must be ''mean'' or ''median''.');
    end

    if nargin < 5 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end

    if nargin < 4 || isempty(filePath)
        filePath = getfield_def(parameters,'reconstructionDir','');
        filePath = fullfile(filePath, 'kcatData');
        if isempty(filePath)
            error('filePath is required or set parameters.reconstructionDir.');
        end
    end

    if nargin < 1 || isempty(DeepLearningModel)
        DeepLearningModel = {'DLKcat','UniKP','CatPred'};
    end
    if isstring(DeepLearningModel) || ischar(DeepLearningModel)
        DeepLearningModel = cellstr(DeepLearningModel);
    end
    DeepLearningModel = unique(strtrim(DeepLearningModel(:)'));
    validTags = {'DLKcat','UniKP','CatPred'};
    if any(~ismember(DeepLearningModel, validTags))
        error('DeepLearningModel must be a subset of {DLKcat, UniKP, CatPred}.');
    end

    % Normalize complex_name list (for robust equality on ReactionName)
    complexNorm = normalize_label_list(complex_name);

    % ------------------------- Load databases -------------------------
    basePath   = fullfile(findECOMAProot, 'scripts','database');
    BRENDA_csv = fullfile(basePath, 'brenda.csv');
    SABIO_csv  = fullfile(basePath, 'sabio.csv');

    T_bren = read_and_prepare_db(BRENDA_csv);
    T_sab  = read_and_prepare_db(SABIO_csv);
    T_db   = [T_bren; T_sab];

    % Provide a stable unique row ID across the concatenated DB
    T_db.row_uid = (1:height(T_db)).';

    % ------------------------- Load predictions ------------------------
    preds = struct();
    for k = 1:numel(DeepLearningModel)
        tag = DeepLearningModel{k};
        pth = find_prediction_path(filePath, tag);
        preds.(tag) = read_and_prepare_pred(pth, tag);
    end

    % ------------------------- Match per model -------------------------
    MATCH = struct();
    for k = 1:numel(DeepLearningModel)
        tag = DeepLearningModel{k};
        MATCH.(tag) = do_match(preds.(tag), T_db, aggMethod, complexNorm);
    end
end


% ============================ Matching core ============================
function MatchTable = do_match(predTbl, dbTbl, aggMethod, complexNorm)
% do_match
%   Run the hierarchical matching and aggregate duplicates by (ProteinID, metaboliteKey).
%   Aggregation uses MEAN or MEDIAN in log10-space as specified by aggMethod.
%   'complexNorm' is a normalized list of complex reaction names for flagging.

    if height(predTbl) == 0 || height(dbTbl) == 0
        MatchTable = table(); return;
    end
    if ~ismember('row_uid', dbTbl.Properties.VariableNames)
        dbTbl.row_uid = (1:height(dbTbl)).';
    end

    % Pre-extract DB columns for speed
    db_uniprot = dbTbl.uniprot;
    db_seq     = dbTbl.sequence;
    db_inchi   = dbTbl.InChIKey;
    db_mnx     = dbTbl.MetaNetXID;

    % Ensure Substrate_norm exists (DB + Pred)
    if ~ismember('Substrate_norm', dbTbl.Properties.VariableNames)
        if ismember('Substrate', dbTbl.Properties.VariableNames)
            dbTbl.Substrate_norm = normalize_substrate_name(dbTbl.Substrate);
        else
            dbTbl.Substrate_norm = strings(height(dbTbl),1);
        end
    end
    db_nameN = dbTbl.Substrate_norm;

    if ~ismember('Substrate_norm', predTbl.Properties.VariableNames)
        if ismember('Substrate', predTbl.Properties.VariableNames)
            predTbl.Substrate_norm = normalize_substrate_name(predTbl.Substrate);
        else
            predTbl.Substrate_norm = strings(height(predTbl),1);
        end
    end

    pred_to_db_predIdx = [];
    pred_to_db_dbIdx   = [];

    % ---------------------- Hierarchical matching ----------------------
    for i = 1:height(predTbl)
        pid = predTbl.ProteinID(i);
        if pid == "" || ismissing(pid), continue; end

        % (1) ProteinID == uniprot
        cand_idx = find(db_uniprot == pid);
        if isempty(cand_idx), continue; end

        % (2) sequence identical IF pred sequence is non-empty
        p_seq = predTbl.sequence(i);
        if ~(p_seq == "" || ismissing(p_seq))
            cand_idx = cand_idx(db_seq(cand_idx) == p_seq);
            if isempty(cand_idx), continue; end
        end

        % (3) Any of {InChIKey, MetaNetXID, Substrate_norm}
        p_inchi = get_field_safe(predTbl, 'InChIKey',       i);
        p_mnx   = get_field_safe(predTbl, 'MetaNetXID',     i);
        p_nameN = get_field_safe(predTbl, 'Substrate_norm', i);

        hit_idx = [];
        for j = cand_idx(:).'
            ok_inchi = (p_inchi ~= "" && p_inchi == db_inchi(j));
            ok_mnx   = (p_mnx   ~= "" && p_mnx   == db_mnx(j));
            ok_name  = (p_nameN ~= "" && p_nameN == db_nameN(j));
            if ok_inchi || ok_mnx || ok_name
                hit_idx(end+1) = j; %#ok<AGROW>
            end
        end

        hit_idx = unique(hit_idx, 'stable'); % avoid triple counting the same DB row
        if isempty(hit_idx), continue; end

        pred_to_db_predIdx = [pred_to_db_predIdx; repmat(i, numel(hit_idx),1)]; %#ok<AGROW>
        pred_to_db_dbIdx   = [pred_to_db_dbIdx;   hit_idx(:)];                   %#ok<AGROW>
    end

    if isempty(pred_to_db_predIdx)
        MatchTable = table(); return;
    end

    % ---------------------- Build raw matched table --------------------
    pred_sub = predTbl(pred_to_db_predIdx, :);
    db_sub   = dbTbl(pred_to_db_dbIdx,  :);

    % Rename DB-side columns to avoid collisions, and carry 'ec'
    db_sub = safe_rename(db_sub, 'InChIKey',       'InChIKey_db');
    db_sub = safe_rename(db_sub, 'MetaNetXID',     'MetaNetXID_db');
    db_sub = safe_rename(db_sub, 'sequence',       'sequence_db');
    db_sub = safe_rename(db_sub, 'value',          'exp_kcat');
    db_sub = safe_rename(db_sub, 'value_log10',    'exp_kcat_log10');
    db_sub = safe_rename(db_sub, 'Substrate',      'Substrate_db');
    db_sub = safe_rename(db_sub, 'Substrate_norm', 'Substrate_norm_db');
    db_sub = safe_rename(db_sub, 'row_uid',        'db_row_uid');
    if ismember('ec', db_sub.Properties.VariableNames)
        db_sub = safe_rename(db_sub, 'ec', 'ec_db');
    else
        db_sub.ec_db = strings(height(db_sub),1);
    end

    % Horizontal concat (ReactionName already in pred_sub if present)
    M = [pred_sub db_sub];

    % ---------------------- Build metaboliteKey (internal) -------------
    H = height(M);
    metKey = strings(H,1);

    hasPredIK  = ismember('InChIKey',     M.Properties.VariableNames) & M.InChIKey~="";
    hasDBIK    = ismember('InChIKey_db',  M.Properties.VariableNames) & M.InChIKey_db~="";
    hasPredMNX = ismember('MetaNetXID',   M.Properties.VariableNames) & M.MetaNetXID~="";
    hasDBMNX   = ismember('MetaNetXID_db',M.Properties.VariableNames) & M.MetaNetXID_db~="";
    hasPredNM  = ismember('Substrate_norm',    M.Properties.VariableNames) & M.Substrate_norm~="";
    hasDBNM    = ismember('Substrate_norm_db', M.Properties.VariableNames) & M.Substrate_norm_db~="";

    % IK priority
    useIK_pred = hasPredIK;
    useIK_db   = ~useIK_pred & hasDBIK;
    metKey(useIK_pred) = "IK:"  + M.InChIKey(useIK_pred);
    metKey(useIK_db)   = "IK:"  + M.InChIKey_db(useIK_db);

    % MNX fallback
    stillEmpty   = metKey == "";
    useMNX_pred  = stillEmpty & hasPredMNX;
    useMNX_db    = stillEmpty & ~hasPredMNX & hasDBMNX;
    metKey(useMNX_pred) = "MNX:" + M.MetaNetXID(useMNX_pred);
    metKey(useMNX_db)   = "MNX:" + M.MetaNetXID_db(useMNX_db);

    % Name fallback
    stillEmpty   = metKey == "";
    useNM_pred   = stillEmpty & hasPredNM;
    useNM_db     = stillEmpty & ~hasPredNM & hasDBNM;
    metKey(useNM_pred) = "NM:" + M.Substrate_norm(useNM_pred);
    metKey(useNM_db)   = "NM:" + M.Substrate_norm_db(useNM_db);

    % ---------------------- Aggregate by group -------------------------
    if ~ismember('predicted_kcat_log10', M.Properties.VariableNames)
        error('do_match:MissingColumn','predicted_kcat_log10 is missing in matched table.');
    end
    if ~ismember('exp_kcat_log10', M.Properties.VariableNames)
        error('do_match:MissingColumn','exp_kcat_log10 is missing in matched table.');
    end

    [G, gProtein, ~] = findgroups(M.ProteinID, metKey);

    x    = M.predicted_kcat_log10;
    y    = M.exp_kcat_log10;
    uids = M.db_row_uid;

    RN_pred = [];
    if ismember('ReactionName', M.Properties.VariableNames)
        RN_pred = M.ReactionName; RN_pred(ismissing(RN_pred)) = "";
        if ~isstring(RN_pred), RN_pred = string(RN_pred); end
    end

    % Original (unnormalized) substrate names for display
    Sub_pred = strings(height(M),1);
    if ismember('Substrate', M.Properties.VariableNames),       Sub_pred = M.Substrate; end
    Sub_db   = strings(height(M),1);
    if ismember('Substrate_db', M.Properties.VariableNames),    Sub_db   = M.Substrate_db; end

    MNX_pred = []; IK_pred = [];
    if ismember('MetaNetXID', M.Properties.VariableNames),     MNX_pred = M.MetaNetXID; end
    if ismember('InChIKey',   M.Properties.VariableNames),     IK_pred  = M.InChIKey;   end

    MNX_db = []; IK_db = [];
    if ismember('MetaNetXID_db', M.Properties.VariableNames),  MNX_db = M.MetaNetXID_db; end
    if ismember('InChIKey_db',   M.Properties.VariableNames),  IK_db  = M.InChIKey_db;   end

    EC_db = M.ec_db;

    nGroups = max(G);
    out_Reaction  = strings(nGroups,1);
    out_Organism  = strings(nGroups,1);
    out_ProteinID = strings(nGroups,1);
    out_ec        = strings(nGroups,1);
    out_Substrate = strings(nGroups,1);
    out_MNX       = strings(nGroups,1);
    out_InChI     = strings(nGroups,1);
    out_n         = zeros(nGroups,1);
    out_x         = NaN(nGroups,1);
    out_y         = NaN(nGroups,1);
    out_isComplex = false(nGroups,1);

    useMedian = strcmpi(aggMethod,'median');
    reducer   = @(v) (useMedian * median(v) + (~useMedian) * mean(v));

    for g = 1:nGroups
        rows = find(G == g);

        xv  = x(rows);
        yv  = y(rows);
        uu  = uids(rows);

        finiteMask = isfinite(xv) & isfinite(yv);
        xv  = xv(finiteMask);
        yv  = yv(finiteMask);
        uu  = uu(finiteMask);

        out_ProteinID(g) = gProtein(g);
        out_n(g)         = numel(unique(uu)); % count UNIQUE DB rows in this group

        % Aggregate log10-values
        if ~isempty(xv), out_x(g) = reducer(xv); end
        if ~isempty(yv), out_y(g) = reducer(yv); end

        % Representative identifiers (prefer pred-side, then DB-side)
        repMNX  = pick_first_nonempty([get_vals(MNX_pred, rows); get_vals(MNX_db, rows)]);
        repIK   = pick_first_nonempty([get_vals(IK_pred,  rows); get_vals(IK_db,  rows)]);
        out_MNX(g)   = repMNX;
        out_InChI(g) = repIK;

        % Substrate for display (prefer original predicted Substrate, else DB)
        repSub = pick_first_nonempty([get_vals(Sub_pred, rows); get_vals(Sub_db, rows)]);
        out_Substrate(g) = repSub;

        % ReactionName from predictions (mode of non-empty)
        if ~isempty(RN_pred)
            out_Reaction(g) = mode_nonempty_string(get_vals(RN_pred, rows));
        end
        
        % Organism from predictions (mode of non-empty)
        if ismember('Organism', M.Properties.VariableNames)
            org_vals = get_vals(M.Organism, rows);
            if ~isstring(org_vals), org_vals = string(org_vals); end
            out_Organism(g) = mode_nonempty_string(org_vals);
        end

        % EC number from DB (mode of non-empty)
        out_ec(g) = mode_nonempty_string(get_vals(EC_db, rows));

        % Complex flag by ReactionName
        rn_norm = normalize_label(out_Reaction(g));
        out_isComplex(g) = any(strcmp(rn_norm, complexNorm));
    end

    % Final aggregated table (NO metaboliteKey; show human-readable Substrate)
    MatchTable = table( ...
        out_Reaction, ...        % ReactionName
        out_Organism, ...        % Organism
        out_ProteinID, ...       % ProteinID
        out_ec, ...              % ec
        out_MNX, ...             % MetaNetXID
        out_InChI, ...           % InChIKey
        out_Substrate, ...       % Substrate
        out_n, ...               % n_group
        out_x, ...               % predicted_kcat_log10
        out_y, ...               % exp_kcat_log10
        out_isComplex, ...       % isComplex
        'VariableNames', {'ReactionName','Organism','ProteinID','ec','MetaNetXID','InChIKey','Substrate', ...
                          'n_group','predicted_kcat_log10','exp_kcat_log10','isComplex'});


    % Linear-space representatives (10.^log10)
    MatchTable.predicted_kcat = 10.^MatchTable.predicted_kcat_log10;
    MatchTable.exp_kcat       = 10.^MatchTable.exp_kcat_log10;
end

% ============================ Utilities ================================
function pth = find_prediction_path(filePath, tag)
    cand = fullfile(filePath, [tag,'.csv']);
    if isfile(cand), pth = cand;
    else, error('BuildKcatMatches:MissingPredFile','File not found: %s', cand);
    end
end

function T = safe_rename(T, oldName, newName)
    if ismember(oldName, T.Properties.VariableNames) && ~strcmp(oldName, newName)
        try
            T = renamevars(T, oldName, newName); % R2020a+
        catch
            T.(newName) = T.(oldName);
            T.(oldName) = [];
        end
    end
end

function val = getfield_def(S, fname, def)
    if ~isstruct(S), val = def; return; end
    if ~isfield(S, fname) || isempty(S.(fname)), val = def; else, val = S.(fname); end
end

function sN = normalize_substrate_name(s)
% Normalize metabolite names for robust equality:
% - lowercase, trim
% - strip trailing bracketed notes
% - unify separators to single space
% - collapse multiple spaces
    if ~isstring(s), s = string(s); end
    sN = lower(strtrim(s));
    sN = regexprep(sN, '\s*(\[[^\]]*\]|\([^\)]*\))\s*$', '', 'once');
    sN = regexprep(sN, '[_\-\,;]+', ' ');
    sN = regexprep(sN, '\s+', ' ');
    sN = strtrim(sN);
end

function lab = normalize_label(s)
% Simple label normalization for reaction names:
% - lowercase
% - trim
% - collapse internal whitespace
    if ~isstring(s), s = string(s); end
    lab = lower(strtrim(s));
    lab = regexprep(lab, '\s+', ' ');
end

function L = normalize_label_list(names)
% Normalize a list of labels for membership checks.
    if nargin == 0 || isempty(names)
        L = strings(0,1); return;
    end
    if ischar(names) || isstring(names), names = cellstr(names); end
    names = string(names(:));
    L = normalize_label(names);
    L = unique(L, 'stable');
end

function v = get_field_safe(T, varName, i)
    if ismember(varName, T.Properties.VariableNames)
        v = T.(varName)(i);
        if ~isstring(v), v = string(v); end
        if ismissing(v), v = ""; end
    else
        v = "";
    end
end

function vals = get_vals(col, idx)
    if isempty(col), vals = strings(numel(idx),1); return; end
    vals = col(idx);
    if ~isstring(vals), vals = string(vals); end
    vals(ismissing(vals)) = "";
end

function out = pick_first_nonempty(s)
    if isempty(s), out = ""; return; end
    if ~isstring(s), s = string(s); end
    for k = 1:numel(s)
        if ~(s(k)=="" || ismissing(s(k))), out = s(k); return; end
    end
    out = "";
end

function out = mode_nonempty_string(s)
% Return the most frequent non-empty string (mode). If all empty -> "".
% Ties are broken by first-seen order.
    if isempty(s), out=""; return; end
    if ~isstring(s), s = string(s); end
    s = s(~ismissing(s) & s~="");
    if isempty(s), out=""; return; end
    [u, ~, idx] = unique(s, 'stable');
    counts = accumarray(idx, 1);
    [~, imax] = max(counts);
    out = u(imax);
end
