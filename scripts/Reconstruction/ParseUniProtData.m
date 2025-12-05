function dbStruct = ParseUniProtData(tsvFilePath)
% ParseUniProtData
%   Parse a UniProt-like TSV and build fast alias index for gene names.
%
% Input:
%   tsvFilePath : path to a TSV with 5 columns:
%                 1) UniProt ID
%                 2) Gene names (may contain multiple aliases, e.g. 'b1518 JW1511')
%                 3) EC numbers
%                 4) Molecular weight (numeric)
%                 5) Sequence
%
% Output (fields):
%   ID            : cellstr, UniProt IDs
%   genes         : cellstr, raw gene-name field as in TSV
%   eccodes       : cellstr, EC numbers
%   MW            : double, molecular weight (NaN if missing)
%   seq           : cellstr, protein sequences
%   geneAliases   : cell-of-cellstr, per-row list of aliases (tokenized, lowercased, unique)
%   aliasMap      : containers.Map(lowercase alias -> row-index vector)
%
% Notes:
%   - aliasMap enables robust matching when your model gene is 'b1518' but the TSV field is 'b1518 JW1511'.
%   - Matching should be done case-insensitively via: q = lower(gene); if isKey(dbStruct.aliasMap,q), idx = dbStruct.aliasMap(q); end

    if ~exist(tsvFilePath, 'file')
        warning('ParseUniProtData:FileNotFound', 'Expected UniProt file not found: %s', tsvFilePath);
        dbStruct = struct('ID',{{}},'genes',{{}},'eccodes',{{}},'MW',[],'seq',{{}}, ...
                          'geneAliases',{{}},'aliasMap',containers.Map('KeyType','char','ValueType','any'));
        return;
    end

    fid = fopen(tsvFilePath, 'r');
    cleaner = onCleanup(@() fclose(fid));

    % Robust textscan: quoted strings, tab-delimited, skip header
    parsedData = textscan(fid, '%q %q %q %q %q', ...
        'Delimiter', '\t', 'HeaderLines', 1, 'ReturnOnError', false, 'EndOfLine', '\n');

    % Populate basics
    dbStruct.ID      = parsedData{1};
    dbStruct.genes   = parsedData{2};
    dbStruct.eccodes = parsedData{3};

    % MW may contain non-numeric tokens; str2double yields NaN which is fine
    dbStruct.MW      = str2double(parsedData{4});
    dbStruct.seq     = parsedData{5};

    n = numel(dbStruct.ID);
    dbStruct.geneAliases = cell(n,1);

    % Build alias map (lower-case keys)
    aliasMap = containers.Map('KeyType','char','ValueType','any');

    for i = 1:n
        raw = '';
        if ~isempty(dbStruct.genes{i})
            raw = dbStruct.genes{i};
        end

        % Split on whitespace and common separators: space, comma, semicolon, slash, pipe
        toks = regexp(raw, '[\s,;/\|]+', 'split');

        % Normalize: trim + lowercase, remove empties
        toks = lower(strtrim(toks));
        toks = toks(~cellfun('isempty', toks));

        % Deduplicate within row
        if ~isempty(toks)
            toks = unique(toks, 'stable');
        end

        dbStruct.geneAliases{i} = toks;

        % Populate aliasMap: alias -> [indices...]
        for t = 1:numel(toks)
            key = toks{t};
            if ~isKey(aliasMap, key)
                aliasMap(key) = i;
            else
                % Append if not already last (avoid many duplicates)
                val = aliasMap(key);
                if isempty(val) || val(end) ~= i
                    aliasMap(key) = [val, i];
                end
            end
        end
    end

    dbStruct.aliasMap = aliasMap;

    % Duplicate UniProt ID check (same as original, but more explicit)
    [uniqueIDs, uniqueIdx] = unique(dbStruct.ID, 'stable');
    if numel(uniqueIDs) < n
        dupIdx = setdiff(1:n, uniqueIdx);
        fprintf('[ParseUniProtData] Found %d duplicate UniProt IDs. Manual inspection may be required.\n', numel(dupIdx));
        disp(dbStruct.ID(dupIdx));
    end
end
