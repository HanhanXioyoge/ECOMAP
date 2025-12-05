function T = read_exp_csv(csvPath)
% READ_EXP_CSV
% Read a CSV into a table, then:
%   1) Rename header "uniprot"  -> "ProteinID"  (case-insensitive)
%   2) Rename header "value"    -> "exp_value"  (case-insensitive)
%   3) Add a new column "exp_value_log10" = log10(exp_value) (<=0 or non-numeric -> NaN)
%
% Notes
% - This function preserves all other columns as-is.
% - It safely handles "value" being stored as text; it converts to double.
% - If both "exp_value" and "value" exist, "exp_value" is preferred as the source.
% - If neither "exp_value" nor "value" exists, the function throws a helpful error.
%
% Example
%   T = read_exp_csv('mydata.csv');

    % --------------------------- Read CSV ---------------------------
    % Preserve original header text as much as MATLAB allows and keep text as string
    T = readtable(csvPath, ...
        'TextType', 'string', ...
        'PreserveVariableNames', true);

    % ------------------------ User-config area ----------------------
    % If you ever want to change the rename rules, edit this block only.
    % Each row is {oldName, newName} matched case-insensitively.
    renamePairs = {
        'uniprot', 'ProteinID';   % "uniprot" -> "ProteinID"
        'value',   'exp_value'    % "value"   -> "exp_value"
    };
    % ----------------------------------------------------------------

    % ------------------------ Apply renames -------------------------
    for k = 1:size(renamePairs,1)
        oldName = renamePairs{k,1};
        newName = renamePairs{k,2};
        T = rename_if_exists_caseinsensitive(T, oldName, newName);
    end

    % ------------------ Identify numeric source ---------------------
    % Prefer "exp_value" if present; else fall back to "value".
    vnamesLower = lower(string(T.Properties.VariableNames));
    hasExp = any(vnamesLower == "exp_value");
    hasVal = any(vnamesLower == "value");

    if ~hasExp && ~hasVal
        error(['Neither "exp_value" nor "value" was found in the table after renaming. ', ...
               'Make sure the CSV has one of these columns (case-insensitive).']);
    end

    sourceName = "exp_value";
    if ~hasExp && hasVal
        sourceName = "value";  % rare case: user’s file already has "value" but not "exp_value"
    end

    % ---------------- Convert to double & compute log10 -------------
    % Make a numeric copy (coerce strings/char to double).
    x = T.(sourceName);
    if isstring(x) || ischar(x) || iscellstr(x) || iscell(x)
        x = str2double(string(x));   % string -> double (non-numeric -> NaN)
    end
    x = double(x); % ensure double

    % If the source was "value", also create/overwrite "exp_value" with the numeric copy
    if sourceName == "value"
        T.exp_value = x;                    % normalized numeric column name
        % (Optionally) drop the old "value" column to avoid ambiguity
        % Comment out next line if you prefer to keep it.
        T = drop_var_caseinsensitive(T, "value");
    else
        % Ensure exp_value is numeric (overwrite if not)
        T.exp_value = x;
    end

    % Compute base-10 log, guarding against non-positive values and NaNs
    logCol = nan(height(T), 1);
    good = isfinite(T.exp_value) & T.exp_value > 0;
    logCol(good) = log10(T.exp_value(good));

    T.exp_value_log10 = logCol;

    % ------------------------- Done ---------------------------------
end

% ====== Helper: case-insensitive rename with collision safety ======
function T = rename_if_exists_caseinsensitive(T, oldName, newName)
    v = T.Properties.VariableNames;
    lowerV = lower(string(v));
    idxOld = find(lowerV == lower(string(oldName)), 1, 'first');
    if isempty(idxOld)
        return; % nothing to do
    end

    % If the newName already exists (case-insensitive), prefer keeping the existing newName
    idxNew = find(lowerV == lower(string(newName)), 1, 'first');

    if isempty(idxNew)
        % Safe rename (works across MATLAB versions)
        v{idxOld} = char(newName);
        T.Properties.VariableNames = v;
    else
        % Both old and new exist; do not blindly rename to avoid duplicate headers.
        % (Optional merge strategy: if newName column is entirely missing while old has values,
        % you could copy non-missing values; kept simple here by doing nothing.)
    end
end

% ====== Helper: drop a variable by name (case-insensitive) ======
function T = drop_var_caseinsensitive(T, varName)
    v = T.Properties.VariableNames;
    lowerV = lower(string(v));
    idx = find(lowerV == lower(string(varName)), 1, 'first');
    if ~isempty(idx)
        T(:, idx) = [];
    end
end
