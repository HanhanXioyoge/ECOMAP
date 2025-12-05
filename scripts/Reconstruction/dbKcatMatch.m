function OUT = dbKcatMatch(DeepLearningModel, filePath, parameters)
    % dbKcatMatch
    % Load DB tables (BRENDA/SABIO), load prediction tables (DLKcat/UniKP/CatPred),
    % perform hierarchical match, deduplicate, and compute correlation stats.
    %
    % The result is packed into:
    % OUT.matches  : Table of matched predictions and experimental data.
    % OUT.stats    : Statistics structure containing correlation data per model.
    % OUT.statsTable : A table summarizing the correlation stats for easy comparison.

    % -------------------- Resolve parameters & paths --------------------
    if nargin < 3 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end

    if nargin < 2 || isempty(filePath)
        filePath = getfield_def(parameters,'dataDir','');
        if isempty(filePath)
            error('filePath is required or set parameters.dataDir.');
        end
    end

    % Normalize model list
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

    % DB file paths
    basePath   = fullfile(findECOMAProot, 'scripts','database');
    BRENDA_csv = fullfile(basePath, 'brenda_merged.csv');
    SABIO_csv  = fullfile(basePath, 'sabio_merged.csv');

    % -------------------- Load databases -----------------
    T_bren = read_and_prepare_db(BRENDA_csv);
    T_sab  = read_and_prepare_db(SABIO_csv);
    T_db   = [T_bren; T_sab];

    % -------------------- Load predictions ----------------
    preds = struct();
    for k = 1:numel(DeepLearningModel)
        tag = DeepLearningModel{k};
        pth = find_prediction_path(filePath, tag);
        preds.(tag) = read_and_prepare_pred(pth, tag);
    end

    % -------------------- Pack raw output -----------------
    OUT = struct();
    OUT.brenda = T_bren;
    OUT.sabio  = T_sab;
    OUT.db     = T_db;
    OUT.preds  = preds;
    OUT.meta   = struct('brenda_csv', BRENDA_csv, ...
                        'sabio_csv',  SABIO_csv,  ...
                        'filePath',   filePath,   ...
                        'models',     {DeepLearningModel});

    % -------------------- Matching & stats ----------------
    OUT.matches = struct();
    OUT.stats   = struct();
    statsData   = []; % Collect stats for all models to form summary table

    for k = 1:numel(DeepLearningModel)
        tag = DeepLearningModel{k};
        thisPred = preds.(tag);

        [tbl, stats] = do_match(thisPred, T_db);

        OUT.matches.(tag) = tbl;
        OUT.stats.(tag)   = stats;

        % Collect stats for comparison
        statsData = [statsData; {tag, stats.corr_r, stats.corr_p, stats.nPoints, stats.MAE_log10}];
    end

    % Create a summary table for easy comparison of models' stats
    statsTable = cell2table(statsData, 'VariableNames', {'Model', 'Correlation_r', 'p_value', 'nPoints', 'MAE_log10'});
    OUT.statsTable = statsTable;
end


% ============================ File finders ==============================
function pth = find_prediction_path(filePath, tag)
    % Prefer exact <tag>.csv; no wildcard fallback here (can add later if needed)
    cand = fullfile(filePath, [tag,'.csv']);
    if isfile(cand)
        pth = cand;
    else
        error('dbKcatMatch:dbMissing', 'The file not found: %s', cand)
    end
end


% ============================ Readers (DB) ==============================
function T = read_and_prepare_db(csvPath)
% read_and_prepare_db
% Read experimental kcat DB CSV (e.g. merged BRENDA / SABIO) and keep only
% the columns we care about for matching + value.
%
% Required columns (must exist in the CSV):
%   uniprot, InChIKey, MetaNetXID, sequence, value
%
% Output columns:
%   uniprot (string)
%   InChIKey (string)
%   MetaNetXID (string)
%   sequence (string)
%   value (double)
%   value_log10 (double, log10(value))

    if ~isfile(csvPath)
        error('dbKcatMatch:dbMissing','DB file not found: %s', csvPath);
    end

    colsAsString = {'uniprot','InChIKey','MetaNetXID','sequence'};
    keepCols     = {'uniprot','InChIKey','MetaNetXID','sequence','value'};

    opts = detectImportOptions(csvPath, 'VariableNamingRule','preserve', ...
                                         'TextType','string');

    % Strong check: these are truly required for downstream matching logic
    missingCols = setdiff(keepCols, opts.VariableNames);
    if ~isempty(missingCols)
        error('read_and_prepare_db:MissingColumns', ...
              'File missing required column(s): %s', strjoin(missingCols, ', '));
    end

    % String columns as string
    presentStringCols = intersect(colsAsString, opts.VariableNames);
    if ~isempty(presentStringCols)
        opts = setvartype(opts, presentStringCols, 'string');
        opts = setvaropts(opts, presentStringCols, 'WhitespaceRule','preserve');
    end

    % Numeric value column
    opts = setvartype(opts, 'value', 'double');

    % Read only the columns we actually need
    opts.SelectedVariableNames = keepCols;
    T = readtable(csvPath, opts);

    % Fill missing strings with ""
    for i = 1:numel(colsAsString)
        v = colsAsString{i};
        if ismember(v, T.Properties.VariableNames)
            if ~isstring(T.(v))
                T.(v) = string(T.(v));
            end
            m = ismissing(T.(v));
            if any(m), T.(v)(m) = ""; end
        end
    end

    % Append log10(value)
    T.value_log10 = log10(T.value);
end


% ============================ Readers (Pred) ============================
function T = read_and_prepare_pred(csvPath, tag)
    % read_and_prepare_pred
    % Read predicted kcat CSV, normalize column names/types, ensure we always
    % have ProteinID / InChIKey / MetaNetXID / sequence columns (string),
    % and predicted_kcat / predicted_kcat_log10 (double).

    if ~isfile(csvPath)
        error('dbKcatMatch:dbMissing','DB file not found: %s', csvPath);
    end

    tag = string(tag);

    % Which column is the numeric prediction?
    switch tag
        case "CatPred"
            predColName  = 'Prediction_(s^(-1))';
        case {"DLKcat","UniKP"}
            predColName  = 'predicted_kcat';
        otherwise
            error('read_and_prepare_pred:UnknownTag','Unknown tag: %s', tag);
    end

    % Columns we *want* to read if present
    wantCols = { ...
        'ProteinID', ...
        'sequence', ...
        'InChIKey', ...
        'MetaNetXID', ...
        predColName ...
    };

    opts = detectImportOptions(csvPath, ...
        'VariableNamingRule','preserve', ...
        'TextType','string');

    % Mandatory columns check: ProteinID + predicted value
    mandatoryCols = {'ProteinID', predColName};
    missingMandatory = setdiff(mandatoryCols, opts.VariableNames);
    if ~isempty(missingMandatory)
        error('read_and_prepare_pred:MissingColumns', ...
              'File missing required column(s): %s', strjoin(missingMandatory, ', '));
    end

    % Keep only what exists in file
    keepCols = intersect(wantCols, opts.VariableNames, 'stable');
    opts.SelectedVariableNames = keepCols;

    % String columns we care about
    stringCols = intersect( ...
        {'ProteinID','sequence','InChIKey','MetaNetXID'}, ...
        keepCols, 'stable');
    if ~isempty(stringCols)
        opts = setvartype(opts, stringCols, 'string');
        opts = setvaropts(opts, stringCols, 'WhitespaceRule','preserve');
    end

    % Prediction column numeric
    opts = setvartype(opts, predColName, 'double');
    opts = setvaropts(opts, predColName, ...
        'TreatAsMissing', {'None','','NA','NaN','nan','N/A'});

    % Read table from file
    T = readtable(csvPath, opts);

    % Drop rows with missing / non-finite prediction
    pv = T.(predColName);
    validMask = ~ismissing(pv) & isfinite(pv);
    T = T(validMask, :);

    % Normalize pred column name to `predicted_kcat`
    unifiedPredVar = 'predicted_kcat';
    if ~strcmp(predColName, unifiedPredVar)
        T = safe_rename(T, predColName, unifiedPredVar);
    end

    % --- Ensure required string columns exist even if missing in file ---
    mustHaveStrCols = {'sequence','InChIKey','MetaNetXID'};
    for c = mustHaveStrCols
        cname = c{1};
        if ~ismember(cname, T.Properties.VariableNames)
            % create empty string column (same height)
            T.(cname) = repmat("", height(T), 1);
        end
    end

    % Fill missing in string columns with ""
    fillCols = {'ProteinID','sequence','InChIKey','MetaNetXID'};
    for i = 1:numel(fillCols)
        v = fillCols{i};
        if ismember(v, T.Properties.VariableNames)
            if ~isstring(T.(v))
                T.(v) = string(T.(v));
            end
            m = ismissing(T.(v));
            if any(m), T.(v)(m) = ""; end
        end
    end

    % Append log10(predicted_kcat)
    if ismember(unifiedPredVar, T.Properties.VariableNames)
        T.([unifiedPredVar '_log10']) = log10(T.(unifiedPredVar));
    else
        % Shouldn't happen, but keep code safe:
        T.predicted_kcat_log10 = NaN(height(T),1);
    end

    % Reorder columns to a nice canonical order
    baseOrder = {'ProteinID','sequence','InChIKey','MetaNetXID', ...
                 'predicted_kcat','predicted_kcat_log10'};
    keepFinal = intersect(baseOrder, T.Properties.VariableNames, 'stable');
    T = T(:, keepFinal);
end


% ============================ Matching helper ==========================
function [MatchTable, MatchStats] = do_match(predTbl, dbTbl)
    % do_match
    % Hierarchical match:
    %   (1) ProteinID == uniprot
    %   (2) sequence identical (only enforced if predTbl.sequence is non-empty)
    %   (3) (InChIKey matches) OR (MetaNetXID matches)
    %
    % Then deduplicate:
    %   For duplicated (ProteinID, MetaNetXID, InChIKey),
    %   keep only the row with minimal |predicted_kcat_log10 - exp_kcat_log10|.
    %
    % Finally compute stats:
    %   - Pearson correlation between predicted_kcat_log10 and exp_kcat_log10
    %   - p-value for the correlation
    %   - nPoints used
    %   - MAE_log10 = mean absolute error in log10 space
    %
    % OUTPUTS
    %   MatchTable : deduplicated best matches (table)
    %   MatchStats : struct with fields corr_r, corr_p, nPoints, MAE_log10

    MatchStats = struct('corr_r', NaN, ...
                        'corr_p', NaN, ...
                        'nPoints', 0, ...
                        'MAE_log10', NaN);

    if height(predTbl) == 0 || height(dbTbl) == 0
        MatchTable = table();
        return;
    end

    % Pre-extract DB columns for speed
    db_uniprot = dbTbl.uniprot;
    db_seq     = dbTbl.sequence;
    db_inchi   = dbTbl.InChIKey;
    db_mnx     = dbTbl.MetaNetXID;

    pred_to_db_predIdx = [];
    pred_to_db_dbIdx   = [];

    for i = 1:height(predTbl)
        pid = predTbl.ProteinID(i);
        if pid == "" || ismissing(pid)
            continue;
        end

        % --- Step 1: ProteinID == uniprot ---
        cand_idx = find(db_uniprot == pid);
        if isempty(cand_idx), continue; end

        % --- Step 2: sequence strict match (if predTbl.sequence is non-empty) ---
        p_seq = predTbl.sequence(i);
        if ~(p_seq == "" || ismissing(p_seq))
            cand_idx = cand_idx(db_seq(cand_idx) == p_seq);
            if isempty(cand_idx), continue; end
        end

        % --- Step 3: metabolite match via InChIKey OR MetaNetXID ---
        p_inchi = predTbl.InChIKey(i);
        p_mnx   = predTbl.MetaNetXID(i);

        hit_idx = [];
        for j = cand_idx(:).'
            ok_inchi = (p_inchi ~= "" && p_inchi == db_inchi(j));
            ok_mnx   = (p_mnx   ~= "" && p_mnx   == db_mnx(j));
            if ok_inchi || ok_mnx
                hit_idx(end+1) = j; %#ok<AGROW>
            end
        end

        if isempty(hit_idx)
            continue;
        end

        % record matches (many-to-many allowed at this stage)
        pred_to_db_predIdx = [pred_to_db_predIdx; repmat(i, numel(hit_idx),1)];
        pred_to_db_dbIdx   = [pred_to_db_dbIdx;   hit_idx(:)];
    end

    % If no hits at all
    if isempty(pred_to_db_predIdx)
        MatchTable = table();
        return;
    end

    % Build raw match table (cartesian pairs that passed the filters)
    pred_sub = predTbl(pred_to_db_predIdx, :);
    db_sub   = dbTbl(pred_to_db_dbIdx,  :);

    % Rename db_sub columns to avoid collisions
    db_sub = safe_rename(db_sub, 'InChIKey',    'InChIKey_db');
    db_sub = safe_rename(db_sub, 'MetaNetXID',  'MetaNetXID_db');
    db_sub = safe_rename(db_sub, 'sequence',    'sequence_db');
    db_sub = safe_rename(db_sub, 'value',       'exp_kcat');
    db_sub = safe_rename(db_sub, 'value_log10', 'exp_kcat_log10');

    % Concat horizontally
    MatchTable = [pred_sub db_sub];

    % ---------------------------------------------------------------------
    % DEDUP STEP:
    % same (ProteinID, MetaNetXID, InChIKey)  → keep row with min |Δlog10|
    % where Δlog10 = predicted_kcat_log10 - exp_kcat_log10
    % ---------------------------------------------------------------------

    % Ensure these columns exist so the grouping key is valid
    if ~ismember('MetaNetXID', MatchTable.Properties.VariableNames)
        MatchTable.MetaNetXID = repmat("", height(MatchTable),1);
    end
    if ~ismember('InChIKey', MatchTable.Properties.VariableNames)
        MatchTable.InChIKey   = repmat("", height(MatchTable),1);
    end

    if ~ismember('predicted_kcat_log10', MatchTable.Properties.VariableNames)
        error('do_match:MissingColumn','predicted_kcat_log10 is missing in MatchTable.');
    end
    if ~ismember('exp_kcat_log10', MatchTable.Properties.VariableNames)
        error('do_match:MissingColumn','exp_kcat_log10 is missing in MatchTable.');
    end

    absDiffAll = abs(MatchTable.predicted_kcat_log10 - MatchTable.exp_kcat_log10);

    % group key = same enzyme & same metabolite definition
    groupKey = MatchTable.ProteinID + "§§" + ...
               MatchTable.MetaNetXID + "§§" + ...
               MatchTable.InChIKey;

    [~, ~, grpIdx] = unique(groupKey, 'stable');

    keepMask = false(height(MatchTable),1);
    for g = 1:max(grpIdx)
        rows = find(grpIdx == g);

        % choose the one with smallest |Δlog10|
        [~, localMinIdx] = min(absDiffAll(rows));
        keepRow = rows(localMinIdx);

        keepMask(keepRow) = true;
    end

    % Deduplicated final table
    MatchTable = MatchTable(keepMask, :);

    % ---------------------------------------------------------------------
    % STATS STEP:
    % Pearson correlation (log10 predicted vs log10 experimental)
    % ---------------------------------------------------------------------
    x = MatchTable.predicted_kcat_log10;
    y = MatchTable.exp_kcat_log10;

    % Only use finite values for stats
    finiteMask = isfinite(x) & isfinite(y);
    x_use = x(finiteMask);
    y_use = y(finiteMask);

    MatchStats.nPoints = numel(x_use);

    if MatchStats.nPoints >= 2
        % corrcoef returns 2x2 matrix: [1 r; r 1]
        C = corrcoef(x_use, y_use, 'Rows','pairwise');
        r = C(1,2);

        % p-value for Pearson r:
        % t = r * sqrt((n-2)/(1-r^2)), df=n-2
        % p = 2 * tcdf(-abs(t), df)
        n  = MatchStats.nPoints;
        t  = r * sqrt((n-2) / max(1e-12, (1-r^2)));
        % MATLAB tcdf needs Statistics Toolbox. If not available, you can skip p.
        % I'll compute p analytically via betainc-free approximation using tcdf,
        % but let's assume tcdf is available. If not, set NaN.
        try
            p = 2 * tcdf(-abs(t), n-2);
        catch
            p = NaN;
        end

        MatchStats.corr_r  = r;
        MatchStats.corr_p  = p;
    else
        MatchStats.corr_r = NaN;
        MatchStats.corr_p = NaN;
    end

    % MAE in log10 space (mean absolute error)
    if MatchStats.nPoints > 0
        MatchStats.MAE_log10 = mean(abs(x_use - y_use));
    else
        MatchStats.MAE_log10 = NaN;
    end
end

% ============================ Small utilities ==========================
function T = safe_rename(T, oldName, newName)
    % Rename table variable oldName -> newName, with fallback for older MATLAB.
    if ismember(oldName, T.Properties.VariableNames)
        if ~strcmp(oldName, newName)
            try
                T = renamevars(T, oldName, newName);  % R2020a+
            catch
                T.(newName) = T.(oldName);
                T.(oldName) = [];
            end
        end
    end
end


function val = getfield_def(S, fname, def)
% getfield_def
% Safe "S.fname or default" helper with warnings.
    if ~isstruct(S)
        warning('dbKcatLoadAll:getfield_def','parameters is not a struct; using default for "%s".', fname);
        val = def; 
        return;
    end

    if ~isfield(S, fname)
        warning('dbKcatLoadAll:getfield_def','parameters.%s not found; using default.', fname);
        val = def; 
        return;
    end

    val = S.(fname);
    if isempty(val)
        warning('dbKcatLoadAll:getfield_def','parameters.%s is empty; using default.', fname);
        val = def;
    end
end
