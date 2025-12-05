function merge_substrate_details(baseCsv, detailsCsv, outCsv)
% MERGE_SUBSTRATE_DETAILS
% Left-join detail columns (all columns after 'substrate' in detailsCsv)
% back into baseCsv by exact match on 'substrate', then place those new
% columns immediately after the 'substrate' column.
%
% Constraints:
%   - Only the 'value' column in baseCsv is numeric; all other columns
%     (including every column from detailsCsv) are treated as strings.
%   - Preserve the original row order of baseCsv.
%   - When writing the output CSV, keep the first column as the ID/index
%     column and leave the first column header blank; keep headers for
%     other columns.
%
% Usage:
%   merge_substrate_details("sabio.csv","sabio_substrate.csv","sabio_merged.csv");

    arguments
        baseCsv    (1,1) string
        detailsCsv (1,1) string
        outCsv     (1,1) string = ""
    end

    if ~isfile(baseCsv);    error('Base file not found: %s', baseCsv);    end
    if ~isfile(detailsCsv); error('Details file not found: %s', detailsCsv); end

    % ---------- Read baseCsv: all string except 'value' kept numeric ----------
    optBase = detectImportOptions(baseCsv, 'VariableNamingRule','preserve');
    allBaseVars = optBase.VariableNames;
    optBase = setvartype(optBase, allBaseVars, 'string');
    if any(strcmp(allBaseVars,'value'))
        optBase = setvartype(optBase, 'value', 'double');
    end
    Tbase = readtable(baseCsv, optBase);

    % ---------- Read detailsCsv: treat ALL columns as string ----------
    optDet = detectImportOptions(detailsCsv, 'VariableNamingRule','preserve');
    optDet = setvartype(optDet, optDet.VariableNames, 'string');
    Tdet   = readtable(detailsCsv, optDet);

    % ---------- Column checks ----------
    if ~ismember('substrate', Tbase.Properties.VariableNames)
        error('Base file is missing the column "substrate" (case-sensitive).');
    end
    if ~ismember('substrate', Tdet.Properties.VariableNames)
        error('Details file is missing the column "substrate" (case-sensitive).');
    end

    % ---------- Trim leading/trailing spaces; keep exact matching ----------
    Tbase.substrate = strtrim(string(Tbase.substrate));
    Tdet.substrate  = strtrim(string(Tdet.substrate));

    % ---------- Preserve row order of baseCsv ----------
    Tbase.rowid_tmp__ = (1:height(Tbase)).';   % temp row-id (must start with a letter)

    % ---------- Select detail columns: all columns after 'substrate' ----------
    vDet = Tdet.Properties.VariableNames;
    idxS = find(strcmp(vDet, 'substrate'), 1, 'first');
    if isempty(idxS); error('Details file: column "substrate" not found.'); end
    if idxS == numel(vDet)
        error('Details file: no columns exist after "substrate" to be merged.');
    end
    detExtraVars = vDet(idxS+1:end);

    % If a substrate appears multiple times in details, keep the first occurrence
    [~, ia] = unique(Tdet.substrate, 'stable');
    Tdet_u  = Tdet(ia, :);

    % Keep only key + detail columns to avoid duplicated columns during join
    Tdet_sel = Tdet_u(:, [{'substrate'}, detExtraVars]);

    % ---------- Left join (base as the left table) ----------
    % 'join' defaults to a left join and keeps only one key column from the left table.
    Tjoined = outerjoin(Tbase, Tdet_sel, ...
        'Keys', 'substrate', ...
        'Type', 'left', ...
        'MergeKeys', true);

    % ---------- Move newly added columns right after 'substrate' ----------
    baseVars  = Tbase.Properties.VariableNames;
    outVars   = Tjoined.Properties.VariableNames;
    addedVars = setdiff(outVars, baseVars, 'stable');  % columns added by join

    idxSub = find(strcmp(outVars, 'substrate'), 1, 'first');
    if isempty(idxSub)
        error('Joined table does not contain the "substrate" column (unexpected).');
    end
    leftPart = outVars(1:idxSub);                      % from col 1 through 'substrate'
    restPart = setdiff(outVars, [leftPart, addedVars], 'stable');
    newOrder = [leftPart, addedVars, restPart];
    Tfinal   = Tjoined(:, newOrder);

    % ---------- Restore original row order; drop temporary row id ----------
    Tfinal = sortrows(Tfinal, 'rowid_tmp__');
    Tfinal.rowid_tmp__ = [];

    % ---------- Write CSV: blank header for the first column; empties as "" ----------
    if outCsv == ""
        [p,b,ext] = fileparts(baseCsv);
        outCsv = fullfile(p, b + "_merged" + ext);
    end
    write_with_blank_first_header(Tfinal, outCsv);

    fprintf('Done: %s\n', outCsv);
end

% ======= Helper: write CSV with blank first header; empty values as "" =======
function write_with_blank_first_header(T, outCsv)
    vars = T.Properties.VariableNames;

    % Convert to cell array for fine-grained control of missing values
    n = height(T); m = width(T);
    C = cell(n, m);

    for j = 1:m
        col = T.(vars{j});

        if isstring(col) || ischar(col)
            s = string(col);
            s(ismissing(s)) = "";               % <missing> -> empty string
            C(:,j) = cellstr(s);

        elseif iscellstr(col) || iscell(col)
            s = strings(n,1);
            for i = 1:n
                if i <= numel(col) && ~isempty(col{i})
                    s(i) = string(col{i});
                else
                    s(i) = "";
                end
            end
            C(:,j) = cellstr(s);

        elseif isnumeric(col) || islogical(col)
            % For numeric columns (e.g., 'value'): NaN -> empty string
            s = strings(n,1);
            if isnumeric(col)
                isn = isnan(col);
                s(~isn) = string(col(~isn));
                s(isn)  = "";
            else
                s = string(col);
            end
            C(:,j) = cellstr(s);

        elseif isdatetime(col) || isduration(col) || iscalendarduration(col)
            s = string(col);
            s(ismissing(s)) = "";
            C(:,j) = cellstr(s);

        else
            s = string(col);
            s(ismissing(s)) = "";
            C(:,j) = cellstr(s);
        end
    end

    % Build header: blank for the first column, keep names for the rest
    header = vars;
    header{1} = "";     % blank first header cell
    ALL = [header; C];

    % writecell will quote fields when needed (commas/quotes/newlines present)
    writecell(ALL, outCsv, 'QuoteStrings', true);
end
