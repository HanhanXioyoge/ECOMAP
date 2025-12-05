function T = read_and_prepare_pred(csvPath, tag)
% READ_AND_PREPARE_PRED
% Read a kcat prediction CSV, normalize column names/types, and ensure
% presence of the canonical columns:
%   ReactionName, Organism, ProteinID, sequence, InChIKey, MetaNetXID,
%   Substrate, Substrate_norm, predicted_kcat, predicted_kcat_log10
%
% Key features:
% - Handles common aliases (case-insensitive), e.g.:
%     * sequence       <- {Protein_sequence, ProteinSequence, AAseq, Seq, ...}
%     * Substrate      <- {Substrate_name, Substrate Name, Ligand, Reactant, ...}
%     * ReactionName   <- {Reaction, rxnName, rxn}
%     * InChIKey       <- {InChI Key, inchi_key, InchiKey}
%     * MetaNetXID     <- {MNX_ID, MNXID, MetaNetX ID}
% - Trims trailing spaces in header names to avoid mismatches like 'sequence '.
% - Drops invalid rows (missing/non-finite/<=0 predicted_kcat).
% - Computes predicted_kcat_log10.
% - Builds Substrate_norm using normalize_substrate_name(Substrate).
%
% Dependencies:
% - normalize_substrate_name (must exist on path)
%
% Inputs
%   csvPath : path to the prediction CSV
%   tag     : "DLKcat" | "UniKP" | "CatPred"
%
% Output
%   T       : normalized table with canonical columns and standard order

    if ~isfile(csvPath)
        error('BuildKcatMatches:predMissing','Prediction file not found: %s', csvPath);
    end

    % -------------------------- Resolve pred column --------------------------
    tag = string(tag);
    switch tag
        case "CatPred"
            predColName  = 'Prediction_(s^(-1))';
        case {"DLKcat","UniKP"}
            predColName  = 'predicted_kcat';
        otherwise
            error('read_and_prepare_pred:UnknownTag','Unknown tag: %s', tag);
    end

    % ------------------------------ Read CSV ---------------------------------
    % We read ALL columns first (do not pre-filter), so that we don't miss
    % columns due to case/space/alias differences. Then we normalize.
    opts = detectImportOptions(csvPath, 'VariableNamingRule','preserve', 'TextType','string');
    % Trim header names to get rid of trailing spaces
    opts.VariableNames = cellstr(strtrim(string(opts.VariableNames)));

    % Make sure the prediction column exists (exact name, case-sensitive here).
    if ~ismember(predColName, opts.VariableNames)
        % Try a case-insensitive fallback for the prediction column name
        actualPred = find_case_insensitive(opts.VariableNames, predColName);
        if isempty(actualPred)
            error('read_and_prepare_pred:MissingColumns','Missing column: %s', predColName);
        else
            predColName = actualPred; % adopt actual name as found
        end
    end

    % Prediction column should be numeric double and treat common missing tokens
    opts = setvartype(opts, predColName, 'double');
    opts = setvaropts(opts, predColName, 'TreatAsMissing', {'None','','NA','NaN','nan','N/A'});

    % Read the full table
    T = readtable(csvPath, opts);
    % Normalize the stored header names (again) post-read
    T.Properties.VariableNames = cellstr(strtrim(string(T.Properties.VariableNames)));

    % ------------------------- Alias normalization ---------------------------
    % Define common aliases (case-insensitive matching).
    aliases.sequence     = {'sequence','Protein_sequence','ProteinSequence','AAseq','Seq'};
    aliases.Substrate    = {'Substrate','Substrate_name','Substrate Name','Ligand','Reactant','substrate_name'};
    aliases.ReactionName = {'ReactionName','Reaction','rxnName','rxn'};
    aliases.InChIKey     = {'InChIKey','InChI Key','inchi_key','InchiKey'};
    aliases.MetaNetXID   = {'MetaNetXID','MNX_ID','MNXID','MetaNetX ID'};

    % If the canonical name is missing but an alias exists, rename the first found alias.
    T = rename_first_match(T, aliases.ReactionName, 'ReactionName');
    T = rename_first_match(T, aliases.sequence,     'sequence');
    T = rename_first_match(T, aliases.InChIKey,     'InChIKey');
    T = rename_first_match(T, aliases.MetaNetXID,   'MetaNetXID');
    T = rename_first_match(T, aliases.Substrate,    'Substrate');

    % ---------------------------- Type cleanup -------------------------------
    % Ensure key string columns exist and are string-typed (fill with "").
    mustString = ["ReactionName","Organism","ProteinID","sequence","InChIKey","MetaNetXID","Substrate"];
    for c = mustString
        if ~ismember(c, T.Properties.VariableNames)
            T.(c) = strings(height(T),1);
        elseif ~isstring(T.(c))
            T.(c) = string(T.(c));
        end
        m = ismissing(T.(c));
        if any(m), T.(c)(m) = ""; end
    end

    % ------------------------ Prediction normalization -----------------------
    % Canonicalize the prediction column name to 'predicted_kcat'
    if ~strcmp(predColName, 'predicted_kcat')
        T = local_rename(T, predColName, 'predicted_kcat');
    end

    % Drop invalid rows: missing / non-finite / <= 0 predicted_kcat
    pv = T.predicted_kcat;
    validMask = ~ismissing(pv) & isfinite(pv) & (pv > 0);
    T = T(validMask, :);

    % Compute log10; safe since pv > 0 after filtering
    T.predicted_kcat_log10 = log10(T.predicted_kcat);

    % -------------------------- Substrate_norm -------------------------------
    % Keep the original Substrate for display; build a normalized variant.
    if ~ismember('Substrate', T.Properties.VariableNames)
        T.Substrate = strings(height(T),1);
    end
    T.Substrate_norm = normalize_substrate_name(T.Substrate);

    % ---------------------------- Column order --------------------------------
    baseOrder = {'ReactionName','Organism','ProteinID','sequence','InChIKey','MetaNetXID', ...
                 'Substrate','Substrate_norm','predicted_kcat','predicted_kcat_log10'};
    keepFinal = intersect(baseOrder, T.Properties.VariableNames, 'stable');
    T = T(:, keepFinal);
end

% ============================== Helpers ======================================

function T = rename_first_match(T, aliasList, canonicalName)
% Rename the first present alias (case-insensitive) to the canonical name,
% only if the canonical name is currently absent.
    if ismember(canonicalName, T.Properties.VariableNames)
        return;
    end
    % Find the first present alias (case-insensitive)
    actual = first_present_case_insensitive(T.Properties.VariableNames, aliasList);
    if ~isempty(actual)
        T = local_rename(T, actual, canonicalName);
    end
end

function nameOut = first_present_case_insensitive(allNames, candidates)
% Return the first existing name in allNames that matches any of candidates
% case-insensitively; otherwise return "".
    allLower = lower(string(allNames));
    for i = 1:numel(candidates)
        idx = find(allLower == lower(string(candidates{i})), 1, 'first');
        if ~isempty(idx)
            nameOut = allNames{idx};
            return;
        end
    end
    nameOut = "";
end

function actual = find_case_insensitive(allNames, target)
% Find the actual header name that case-insensitively matches 'target'.
    idx = find(lower(string(allNames)) == lower(string(target)), 1, 'first');
    if isempty(idx)
        actual = "";
    else
        actual = allNames{idx};
    end
end

function T = local_rename(T, oldName, newName)
% Minimal, version-agnostic variable rename (no toolbox dependency).
    if ~ismember(oldName, T.Properties.VariableNames) || strcmp(oldName, newName)
        return;
    end
    vn = T.Properties.VariableNames;
    vn{strcmp(vn, oldName)} = newName;
    T.Properties.VariableNames = vn;
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