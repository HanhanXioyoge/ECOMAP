function T_all = collapse_match_to_single_table(OUT, varargin)
% COLLAPSE_MATCH_TO_SINGLE_TABLE
%   Collapse OUT.(modelName) tables into ONE table by inner-joining on
%   {'ProteinID','InChIKey'} and renaming the ONLY differing columns
%   'predicted_kcat_log10' & 'predicted_kcat' to suffixed names per model.
%
%   T_all = collapse_match_to_single_table(OUT)
%   T_all = collapse_match_to_single_table(OUT, 'Models', {'DLKcat','UniKP'})
%   T_all = collapse_match_to_single_table(OUT, 'NormalizeInChIKey', true, 'CaseSensitive', true)
%
% INPUT
%   OUT : struct with fields like OUT.DLKcat/OUT.UniKP/...; each field is a table.
%
% NAME-VALUE OPTIONS
%   'Models'            : cellstr of model names to include (default: {'DLKcat','UniKP','CatPred'})
%   'KeyCols'           : cellstr, length 2, default: {'ProteinID','InChIKey'}
%   'PredCols'          : cellstr, exactly {'predicted_kcat_log10','predicted_kcat'} by default
%   'NormalizeInChIKey' : logical, strip 'InChIKey=' prefix and UPPERCASE (default: true)
%   'CaseSensitive'     : logical, key comparison case-sensitive? (default: true)
%
% OUTPUT
%   T_all : single table containing keys and model-suffixed prediction columns.
%
% NOTES
%   - Keeps only rows **present in ALL specified models** (inner join).
%   - Within each model, duplicates by key are collapsed (keep first).
%   - Only the key columns and the two prediction columns are kept.
%   - Prediction columns are renamed to: predicted_kcat_log10_<MODEL>, predicted_kcat_<MODEL>.
%

    % ------------ parse inputs ------------
    ip = inputParser;
    ip.addParameter('Models', {'DLKcat','UniKP','CatPred'}, @(c)iscell(c)||isstring(c));
    ip.addParameter('KeyCols', {'ProteinID','InChIKey', 'Organism', 'ec', 'isComplex', 'exp_kcat_log10', 'exp_kcat'}, @(c)iscell(c)&&numel(c)==7);
    ip.addParameter('PredCols', {'predicted_kcat_log10','predicted_kcat'}, @(c)iscell(c)&&numel(c)==2);
    ip.addParameter('NormalizeInChIKey', true, @(x)islogical(x)&&isscalar(x));
    ip.addParameter('CaseSensitive', true, @(x)islogical(x)&&isscalar(x));
    ip.parse(varargin{:});

    models            = cellstr(string(ip.Results.Models));
    keyCols           = ip.Results.KeyCols;
    predCols          = ip.Results.PredCols;
    normalizeInChIKey = ip.Results.NormalizeInChIKey;
    caseSensitive     = ip.Results.CaseSensitive;

    % ------------ sanity checks ------------
    for k = 1:numel(models)
        mdl = models{k};
        if ~isfield(OUT, mdl) || ~istable(OUT.(mdl))
            error('collapse_match_to_single_table:MissingModelTable', ...
                  'OUT.%s not found or not a table.', mdl);
        end
        T = OUT.(mdl);
        need = [keyCols, predCols];
        miss = setdiff(need, T.Properties.VariableNames);
        if ~isempty(miss)
            error('collapse_match_to_single_table:MissingColumns', ...
                 'OUT.%s missing columns: %s', mdl, strjoin(miss, ', '));
        end
    end

    % ------------ build per-model trimmed tables ------------
    % We keep:
    %   - original key columns from the *first* model for output
    %   - normalized join keys K1/K2 for robust matching
    %   - two prediction columns renamed with suffix _<MODEL>
    base = models{1};
    T_all = prepare_one(OUT.(base), base, keyCols, predCols, normalizeInChIKey, caseSensitive, true);

    for k = 2:numel(models)
        mdl = models{k};
        Tk  = prepare_one(OUT.(mdl), mdl, keyCols, predCols, normalizeInChIKey, caseSensitive, false);

        % inner join on normalized keys; only bring over suffixed pred columns from right
        rightVars = setdiff(Tk.Properties.VariableNames, {'K1','K2', keyCols{:}});
        T_all = innerjoin(T_all, Tk, 'Keys', {'K1','K2'}, ...
                          'LeftVariables', T_all.Properties.VariableNames, ...
                          'RightVariables', rightVars);
    end

    % Drop normalized join keys, keep original key columns once
    T_all = removevars(T_all, intersect(T_all.Properties.VariableNames, {'K1','K2'}));
end

% ======================== helpers ========================

function Tprep = prepare_one(T, modelName, keyCols, predCols, normInChI, caseSensitive, keepOriginalKeys)
% Trim to keys + preds, normalize keys, dedup by key, and rename preds with suffix.

    % Make safe string columns for keys
    PID = to_str_col(T.(keyCols{1}));
    IK  = to_str_col(T.(keyCols{2}));

    PID = strip(PID);
    IK  = strip(IK);

    if normInChI
        IK = regexprep(IK, '^InChIKey\s*=\s*', '', 'ignorecase');
        IK = upper(IK);
    end
    if ~caseSensitive
        PID = lower(PID);
        IK  = lower(IK);
    end

    % Build normalized join keys
    K1 = PID;
    K2 = IK;

    % Deduplicate by normalized key (keep first)
    K = K1 + "|" + K2;
    [~, firstIdx] = unique(K, 'stable');

    % Slice variables to keep
    varsToKeep = [keyCols, predCols];
    Tsub = T(firstIdx, varsToKeep);

    % Attach normalized keys
    Tsub.K1 = K1(firstIdx);
    Tsub.K2 = K2(firstIdx);

    % Optionally drop the original visible key columns for non-base models
    if ~keepOriginalKeys
        Tsub = removevars(Tsub, keyCols);
    end

    % Rename the two prediction columns with model suffix
    suf = ['_' char(modelName)];
    Tsub.Properties.VariableNames( ...
        ismember(Tsub.Properties.VariableNames, predCols)) = ...
        strcat(predCols, suf);

    Tprep = Tsub;
end

function s = to_str_col(col)
% Robustly convert a table column to string array.
    if isstring(col)
        s = col;
    elseif iscellstr(col)
        s = string(col);
    elseif iscell(col)
        s = strings(size(col));
        for i = 1:numel(col)
            if isstring(col{i})
                s(i) = col{i};
            elseif ischar(col{i})
                s(i) = string(col{i});
            elseif isnumeric(col{i}) || islogical(col{i})
                s(i) = string(col{i});
            else
                s(i) = "";
            end
        end
    elseif isnumeric(col) || islogical(col)
        s = string(col);
    else
        s = string(col);
    end
end
