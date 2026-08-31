function candidates = fetchCrossSpeciesEnzymes(seeds, opts)
%FETCHCROSSSPECIESENZYMES Retrieve cross-species enzyme candidates from UniProt.
%
%   For each seed in the input table, fetches UniProt proteins matching
%   the cleaned reaction name (or EC number as fallback), and applies
%   quality filters and per-organism deduplication.
%
%   candidates = fetchCrossSpeciesEnzymes(seeds)
%   candidates = fetchCrossSpeciesEnzymes(seeds, opts)
%
%   Inputs:
%     seeds - table from buildOkoReactionSeeds with columns:
%       ReactionName, OriginalProteinID, ECNumber, ...
%     opts  - struct with fields:
%       MaxHomologs        (default 100)
%       MinLength          (default 30 aa)
%       MaxLength          (default 2000 aa)
%       ExcludeFragment    (default true)
%       RequireReviewed    (default false)
%       MaxRetrievedRecords (default 5000)
%       CacheDir           (default '')
%       RateLimit          (default 0.2 s/request)
%       CheckpointCsv      durable candidatesAll CSV path (default '')
%       CheckpointMat      resume state MAT path (default '')
%       CheckpointEvery    persist after this many completed seeds (25)
%       Resume             reuse a compatible CheckpointMat (default true)
%       CheckpointFn       optional @(candidatesAll, processedKeys) callback
%       MockServer         (default []) — optional mockUniProtServer
%                          instance for testing (skip HTTP)
%
%   Output:
%     candidates - table with columns:
%       ReactionName       - copied from input seed (verbatim)
%       OriginalProteinID  - copied from input seed
%       accession          - homolog UniProt
%       organism           - homolog scientific name
%       taxID              - homolog NCBI tax ID
%       proteinName        - homolog protein name
%       ecNumbers          - homolog EC annotation
%       sequence           - homolog sequence
%       length             - homolog sequence length
%       fragment           - UniProt fragment flag
%       reviewed           - UniProt reviewed flag
%       dateSequenceModified - UniProt date_sequence_modified

    if nargin < 2
        opts = struct();
    end
    opts = applyDefaults(opts);

    % Per-seed fetch + filter + dedup
    candidatesAll = table();
    seedKeys = makeSeedKeys(seeds);
    processedKeys = strings(0,1);
    if opts.Resume && ~isempty(opts.CheckpointMat) && exist(opts.CheckpointMat, 'file')
        try
            state = load(opts.CheckpointMat, 'candidatesAll', 'processedKeys', 'seedKeys');
            if isfield(state, 'seedKeys') && isequal(string(state.seedKeys), seedKeys)
                if isfield(state, 'candidatesAll'), candidatesAll = state.candidatesAll; end
                if isfield(state, 'processedKeys'), processedKeys = string(state.processedKeys); end
                fprintf('[fetchCrossSpeciesEnzymes] Resuming: %d/%d seeds complete\n', ...
                    numel(processedKeys), height(seeds));
            else
                warning('fetchCrossSpeciesEnzymes:checkpointMismatch', ...
                    'Ignoring candidate checkpoint because the seed set changed.');
            end
        catch err
            warning('fetchCrossSpeciesEnzymes:checkpointReadFailed', ...
                'Ignoring unreadable candidate checkpoint: %s', err.message);
        end
    end
    % A genome-scale model commonly contains thousands of reaction arms but
    % far fewer unique (enzyme-name, EC) searches.  Reuse lookups within one
    % run so isozymes/repeated substrates do not hammer UniProt repeatedly.
    lookupCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:height(seeds)
        seedKey = seedKeys(i);
        if any(processedKeys == seedKey), continue; end
        cleanName = cleanEnzymeName(seeds.ReactionName{i});
        ecNumber = seeds.ECNumber{i};
        if strcmp(ecNumber, ''), ecNumber = ''; end

        perOpts = opts;
        perOpts.TargetSeq = '';

        lookupKey = lower(strtrim([cleanName '||' ecNumber]));
        if isKey(lookupCache, lookupKey)
            cand = lookupCache(lookupKey);
        else
            try
                cand = fetchOneEnzyme(cleanName, ecNumber, perOpts);
            catch err
                location = '';
                if ~isempty(err.stack)
                    location = sprintf(' (%s:%d)', err.stack(1).name, err.stack(1).line);
                end
                warning('fetchCrossSpeciesEnzymes:fetchFailed', ...
                        'Failed to fetch %s%s: %s', cleanName, location, err.message);
                cand = table();
                % Do not mark this seed complete. A later invocation will
                % retry it while retaining all successfully processed seeds.
                maybeCheckpoint(false);
                continue;
            end
            lookupCache(lookupKey) = cand;
        end

        if ~isempty(cand)
            % Stamp seed identity onto each homolog row
            nRows = height(cand);
            seedRow = seeds(i, :);
            cand.ReactionName      = repmat({seedRow.ReactionName{1}}, nRows, 1);
            cand.OriginalProteinID = repmat({seedRow.OriginalProteinID{1}}, nRows, 1);
            candidatesAll = [candidatesAll; cand]; %#ok<AGROW>
        end
        processedKeys(end+1,1) = seedKey; %#ok<AGROW>
        maybeCheckpoint(i == height(seeds));
    end

    % Covers an all-resumed run and a final failed seed.
    persistCheckpoint();
    incomplete = ~ismember(seedKeys, processedKeys);
    if any(incomplete)
        error('fetchCrossSpeciesEnzymes:IncompleteCheckpoint', ...
            ['UniProt retrieval stopped with %d unfinished seeds. ' ...
             'All completed candidates were checkpointed; rerun the same ' ...
             'command to continue.'], sum(incomplete));
    end

    candidates = candidatesAll;

    function maybeCheckpoint(force)
        if force || (opts.CheckpointEvery > 0 && ...
                mod(numel(processedKeys), opts.CheckpointEvery) == 0)
            persistCheckpoint();
        end
    end

    function persistCheckpoint()
        if ~isempty(opts.CheckpointCsv) && ~isempty(candidatesAll)
            atomicWriteTable(candidatesAll, opts.CheckpointCsv);
        end
        if ~isempty(opts.CheckpointMat)
            atomicSaveState(opts.CheckpointMat, candidatesAll, processedKeys, seedKeys);
        end
        if ~isempty(opts.CheckpointFn)
            opts.CheckpointFn(candidatesAll, processedKeys);
        end
    end
end

function opts = applyDefaults(opts)
    if ~isfield(opts, 'MaxHomologs'),         opts.MaxHomologs = 100; end
    if ~isfield(opts, 'MinLength'),           opts.MinLength = 30; end
    if ~isfield(opts, 'MaxLength'),           opts.MaxLength = 2000; end
    if ~isfield(opts, 'ExcludeFragment'),     opts.ExcludeFragment = true; end
    if ~isfield(opts, 'RequireReviewed'),     opts.RequireReviewed = false; end
    if ~isfield(opts, 'MaxRetrievedRecords'), opts.MaxRetrievedRecords = 5000; end
    if ~isfield(opts, 'CacheDir'),            opts.CacheDir = ''; end
    if ~isfield(opts, 'RateLimit'),           opts.RateLimit = 0.2; end
    if ~isfield(opts, 'MockServer'),          opts.MockServer = []; end
    if ~isfield(opts, 'CheckpointCsv'),       opts.CheckpointCsv = ''; end
    if ~isfield(opts, 'CheckpointMat'),       opts.CheckpointMat = ''; end
    if ~isfield(opts, 'CheckpointEvery'),     opts.CheckpointEvery = 25; end
    if ~isfield(opts, 'Resume'),              opts.Resume = true; end
    if ~isfield(opts, 'CheckpointFn'),        opts.CheckpointFn = []; end
end

function keys = makeSeedKeys(seeds)
    keys = strings(height(seeds),1);
    for i = 1:height(seeds)
        keys(i) = string(seeds.ReactionName{i}) + "||" + ...
                  string(seeds.OriginalProteinID{i}) + "||" + ...
                  string(seeds.ECNumber{i});
    end
end

function cand = fetchOneEnzyme(cleanName, ecNumber, opts)
% Fetch one enzyme's cross-species candidates.
%
% Uses UniProt REST cursor pagination, deterministic sort, per-organism
% dedup, sequence dedup, and MaxHomologs cap. Implementation mirrors
% the original fetchUniProtHomologs behavior (per spec §4.2).
    cand = table();
    if isempty(opts.MockServer)
        rawRows = httpFetch(cleanName, ecNumber, opts);
        % Protein names in ecModels are not always UniProt canonical names.
        % Fall back to the exact EC query before declaring the seed uncovered.
        if isempty(rawRows) && ~isempty(ecNumber)
            rawRows = httpFetch('', ecNumber, opts);
        end
    else
        rawRows = opts.MockServer.fetch(cleanName, ecNumber);
    end
    if isempty(rawRows), return; end

    rows = filterBySequence(rawRows);
    rows = filterByLength(rows, opts.MinLength, opts.MaxLength);
    if opts.ExcludeFragment
        rows = filterByFragment(rows);
    end
    if opts.RequireReviewed
        rows = filterByReviewed(rows);
    end
    if ~isempty(ecNumber)
        rows = filterByEC(rows, ecNumber);
    end
    rows = deterministicSort(rows, cleanName);
    rows = dedupPerOrganism(rows);
    rows = dedupBySequence(rows);
    if height(rows) > opts.MaxHomologs
        rows = rows(1:opts.MaxHomologs, :);
    end

    cand = rowsToTable(rows);
end

function rawRows = httpFetch(cleanName, ecNumber, opts)
% UniProt REST cursor pagination via Link: rel="next"
    if isempty(cleanName)
        query = sprintf('ec:%s', ecNumber);
    else
        query = sprintf('protein_name:"%s"', cleanName);
    end
    if ~isempty(cleanName) && ~isempty(ecNumber)
        query = sprintf('%s AND ec:%s', query, ecNumber);
    end

    fields = 'accession,organism_name,organism_id,protein_name,ec,sequence,length,fragment,reviewed,date_sequence_modified';
    cacheKey = computeCacheKey(query, fields);

    % Cache check
    cachePath = fullfile(opts.CacheDir, [cacheKey '.tsv']);
    if ~isempty(opts.CacheDir) && exist(cachePath, 'file')
        rawRows = loadFromCache(cachePath);
        if ~isempty(rawRows), return; end
    end

    baseURL = 'https://rest.uniprot.org/uniprotkb/search';
    url = sprintf('%s?query=%s&fields=%s&format=tsv&size=500', ...
                  baseURL, urlencode(query), urlencode(fields));

    rawRows = table();
    totalFetched = 0;
    pageCount = 0;
    maxPages = ceil(opts.MaxRetrievedRecords / 500);

    while ~isempty(url) && pageCount < maxPages && totalFetched < opts.MaxRetrievedRecords
        pageCount = pageCount + 1;
        pause(opts.RateLimit);
        [respBody, linkHeader] = httpGetWithRetry(url, opts);
        if isempty(respBody), break; end

        pageTable = parseTsvPage(respBody);
        rawRows = [rawRows; pageTable];
        totalFetched = height(rawRows);
        url = extractNextLink(linkHeader);
    end

    if ~isempty(opts.CacheDir)
        saveToCache(cachePath, rawRows);
    end
end

function col = findCol(T, candidates)
%FINDCOL Return the first existing column name from the candidates list,
%   or '' if none exist. readtable may rename columns that conflict with
%   builtins or have duplicates; this helper provides fallback aliases.
    col = '';
    for i = 1:numel(candidates)
        if ismember(candidates{i}, T.Properties.VariableNames)
            col = candidates{i};
            return;
        end
    end
end

function v = colOrDefault(T, candidates, default)
%COLORDEFAULT Read column from T by trying candidates in order; if none
%   exist, return default.
    col = findCol(T, candidates);
    if isempty(col)
        v = default;
        return;
    end
    v = T.(col);
end

function [body, linkHeader] = httpGetWithRetry(url, ~)
    maxAttempts = 3;
    backoff = 1;
    body = '';
    linkHeader = '';
    lastError = [];
    for attempt = 1:maxAttempts
        try
            try
                noCompression = matlab.net.http.HeaderField('Accept-Encoding', 'identity');
                rm = matlab.net.http.RequestMessage('GET', noCompression);
                httpOptions = matlab.net.http.HTTPOptions('ConnectTimeout', 60);
                resp = rm.send(url, httpOptions);
                if double(resp.StatusCode) < 200 || double(resp.StatusCode) >= 300
                    error('fetchCrossSpeciesEnzymes:HttpStatus', ...
                        'UniProt returned HTTP %d.', double(resp.StatusCode));
                end
                body = decodeHttpBody(resp.Body.Data);
                linkFields = resp.getFields('Link');
                if ~isempty(linkFields)
                    linkHeader = char(string(linkFields(1).Value));
                end
                return;
            catch
                body = webread(url);
                return;
            end
        catch err
            lastError = err;
            if attempt == maxAttempts
                rethrow(lastError);
            end
            pause(backoff);
            backoff = backoff * 2;
        end
    end
end

function body = decodeHttpBody(data)
%DECODEHTTPBODY Decode UniProt text responses, including gzip payloads.
    if ischar(data), body = data; return; end
    if isstring(data), body = char(data); return; end
    bytes = uint8(data(:));
    if numel(bytes) >= 2 && bytes(1) == hex2dec('1F') && bytes(2) == hex2dec('8B')
        gzPath = [tempname '.gz'];
        fid = fopen(gzPath, 'wb');
        fwrite(fid, bytes, 'uint8');
        fclose(fid);
        cleanupGz = onCleanup(@()deleteIfPresent(gzPath)); %#ok<NASGU>
        outDir = tempname; mkdir(outDir);
        cleanupDir = onCleanup(@()removeDirIfPresent(outDir)); %#ok<NASGU>
        files = gunzip(gzPath, outDir);
        if iscell(files), decodedPath = files{1}; else, decodedPath = char(files(1)); end
        fid = fopen(decodedPath, 'rb');
        raw = fread(fid, '*uint8')';
        fclose(fid);
        body = native2unicode(raw, 'UTF-8');
    else
        body = native2unicode(bytes', 'UTF-8');
    end
end

function deleteIfPresent(path)
    if exist(path, 'file'), delete(path); end
end

function removeDirIfPresent(path)
    if exist(path, 'dir'), rmdir(path, 's'); end
end

function nextURL = extractNextLink(linkHeader)
    nextURL = '';
    if isempty(linkHeader), return; end
    tokens = regexp(linkHeader, '<([^>]+)>;\s*rel="next"', 'tokens', 'once');
    if ~isempty(tokens)
        nextURL = tokens{1};
    end
end

function rows = filterBySequence(rawRows)
    seqCol = findCol(rawRows, {'sequence', 'sequence_1'});
    if isempty(seqCol)
        rows = rawRows;
        return;
    end
    seqs = rawRows.(seqCol);
    mask = strlength(strtrim(string(seqs))) > 0;
    rows = rawRows(mask, :);
end

function rows = filterByLength(rawRows, mn, mx)
    lenCol = findCol(rawRows, {'length'});
    if isempty(lenCol)
        rows = rawRows;
        return;
    end
    rawLengths = rawRows.(lenCol);
    if isnumeric(rawLengths), lens = double(rawLengths); else, lens = str2double(string(rawLengths)); end
    mask = lens >= mn & lens <= mx;
    rows = rawRows(mask, :);
end

function rows = filterByFragment(rawRows)
    fragCol = findCol(rawRows, {'fragment'});
    if isempty(fragCol)
        rows = rawRows;
        return;
    end
    frag = rawRows.(fragCol);
    if islogical(frag)
        mask = ~frag;
    elseif isnumeric(frag)
        % An entirely empty UniProt fragment column is inferred by
        % readtable as NaN doubles. Empty means "not annotated fragment".
        mask = isnan(frag) | frag == 0;
    else
        values = textVector(frag);
        mask = ~(strcmpi(strtrim(values), "true") | ...
                 strcmpi(strtrim(values), "fragment") | ...
                 strcmpi(strtrim(values), "1"));
    end
    rows = rawRows(mask, :);
end

function rows = filterByReviewed(rawRows)
    revCol = findCol(rawRows, {'reviewed'});
    if isempty(revCol), rows = rawRows([],:); return; end
    rev = rawRows.(revCol);
    if islogical(rev)
        mask = rev;
    else
        values = textVector(rev);
        mask = strcmpi(strtrim(values), 'true') | ...
               strcmpi(strtrim(values), 'reviewed');
    end
    rows = rawRows(mask, :);
end

function rows = filterByEC(rawRows, ecNumber)
    ecCol = findCol(rawRows, {'ec'});
    if isempty(ecCol)
        rows = rawRows;
        return;
    end
    ec = rawRows.(ecCol);
    ecText = textVector(ec);
    mask = false(height(rawRows), 1);
    for i = 1:height(rawRows)
        v = char(ecText(i));
        % Exact match against EC list (comma/semicolon-separated)
        if ~isempty(v)
            tokens = regexp(v, '[;,]\s*', 'split');
            for t = 1:numel(tokens)
                if strcmp(strtrim(tokens{t}), ecNumber)
                    mask(i) = true;
                    break;
                end
            end
        end
    end
    rows = rawRows(mask, :);
end

function rows = deterministicSort(rawRows, cleanName)
% Sort by (reviewed desc, exact-name desc, length-closeness desc, accession asc).
    targetLen = 0;
    seqCol = findCol(rawRows, {'sequence', 'sequence_1'});
    if ~isempty(seqCol)
        seqs = string(rawRows.(seqCol));
        lens = strlength(seqs);
        if ~isempty(lens), targetLen = median(lens); end
    end
    revCol = findCol(rawRows, {'reviewed'});
    if isempty(revCol)
        reviewedScore = zeros(height(rawRows), 1);
    elseif islogical(rawRows.(revCol))
        reviewedScore = double(rawRows.(revCol));
    else
        values = textVector(rawRows.(revCol));
        reviewedScore = double(strcmpi(strtrim(values), "true") | ...
                               strcmpi(strtrim(values), "reviewed"));
    end
    pnCol = findCol(rawRows, {'protein_name'});
    if isempty(pnCol)
        exactScore = zeros(height(rawRows), 1);
    else
        exactScore = double(strcmpi(rawRows.(pnCol), cleanName));
    end
    if isempty(seqCol)
        lens = zeros(height(rawRows), 1);
    else
        lens = double(strlength(string(rawRows.(seqCol))));
    end
    lenScore = -abs(lens - targetLen);

    accCol = findCol(rawRows, {'accession'});
    if isempty(accCol)
        accScore = zeros(height(rawRows), 1);
    else
        accs = rawRows.(accCol);
        [~, ~, accScore] = unique(accs);
    end

    key = [-reviewedScore, -exactScore, -lenScore, accScore];
    [~, order] = sortrows(key);
    rows = rawRows(order, :);
end

function rows = dedupPerOrganism(rawRows)
% Keep one row per organism. After deterministicSort, the reviewed-first
% row for each organism is at the top.
    orgCol = findCol(rawRows, {'organism_name'});
    if isempty(orgCol)
        rows = rawRows;
        return;
    end
    orgs = rawRows.(orgCol);
    [~, ia] = unique(orgs, 'stable');
    rows = rawRows(ia, :);
end

function rows = dedupBySequence(rawRows)
% Drop rows with duplicate sequences (exact string match).
    seqCol = findCol(rawRows, {'sequence', 'sequence_1'});
    if isempty(seqCol)
        rows = rawRows;
        return;
    end
    seqs = string(rawRows.(seqCol));
    [~, ia] = unique(seqs, 'stable');
    rows = rawRows(ia, :);
end

function T = rowsToTable(rows)
%ROWSTOTABLE Convert UniProt TSV rows to the canonical candidates table.
%
%   Uses findCol / colOrDefault for defensive column access. readtable
%   may rename columns that conflict with builtins or have duplicates;
%   the helpers here try a list of aliases before falling back to defaults.
    T = table();
    T.accession = cellstr(textVector(colOrDefault(rows, {'accession'}, repmat({''}, height(rows), 1))));
    T.organism = cellstr(textVector(colOrDefault(rows, {'organism_name'}, repmat({''}, height(rows), 1))));
    T.taxID = numericVector(colOrDefault(rows, {'organism_id'}, zeros(height(rows), 1)));
    T.proteinName = cellstr(textVector(colOrDefault(rows, {'protein_name'}, repmat({''}, height(rows), 1))));
    T.ecNumbers = cellstr(textVector(colOrDefault(rows, {'ec'}, repmat({''}, height(rows), 1))));
    T.sequence = cellstr(textVector(colOrDefault(rows, {'sequence', 'sequence_1'}, repmat({''}, height(rows), 1))));
    T.length = numericVector(colOrDefault(rows, {'length'}, zeros(height(rows), 1)));

    fragCol = findCol(rows, {'fragment'});
    if isempty(fragCol)
        T.fragment = false(height(rows), 1);
    else
        frag = rows.(fragCol);
        if islogical(frag)
            T.fragment = frag;
        elseif isnumeric(frag)
            T.fragment = ~(isnan(frag) | frag == 0);
        else
            values = textVector(frag);
            T.fragment = strcmpi(strtrim(values), 'true') | ...
                         strcmpi(strtrim(values), 'fragment') | ...
                         strcmpi(strtrim(values), '1');
        end
    end

    revCol = findCol(rows, {'reviewed'});
    if isempty(revCol)
        T.reviewed = false(height(rows), 1);
    else
        rev = rows.(revCol);
        if islogical(rev)
            T.reviewed = rev;
        else
            values = textVector(rev);
            T.reviewed = strcmpi(strtrim(values), 'true') | ...
                         strcmpi(strtrim(values), 'reviewed');
        end
    end

    T.dateSequenceModified = cellstr(textVector(colOrDefault(rows, ...
        {'date_sequence_modified', 'date_modified'}, ...
        repmat({''}, height(rows), 1))));
end

function values = textVector(values)
%TEXTVECTOR Convert heterogeneous external text columns without allowing
% <missing>, undefined categoricals or NaT values to reach char/cellstr.
    if ischar(values)
        % A table variable can be an N-by-M character matrix, with one text
        % value per row. Preserve those rows instead of flattening chars.
        values = string(cellstr(values));
    else
        values = string(values(:));
    end
    values(ismissing(values)) = "";
end

function values = numericVector(values)
    if isnumeric(values)
        values = double(values(:));
    else
        values = str2double(string(values(:)));
    end
end

function T = parseTsvPage(body)
%PARSEPAGE Parse a UniProt TSV response body into a table.
%
%   readtable interprets a string argument as a filename, so we write
%   the body to a temp file first and parse from there. Use
%   'VariableNamingRule','modify' so readtable silently renames any
%   duplicate or invalid column names instead of emitting warnings.
    % Replace the external header before readtable sees it.  Some UniProt
    % display headers exceed namelengthmax and MATLAB warns even when
    % VariableNamingRule='preserve'.
    firstBreak = regexp(body, '\r?\n', 'once');
    if isempty(firstBreak), T = table(); return; end
    dataBody = body(firstBreak+1:end);
    canonical = {'accession','organism_name','organism_id','protein_name', ...
                 'ec','sequence','length','fragment','reviewed', ...
                 'date_sequence_modified'};
    normalizedBody = [strjoin(canonical, sprintf('\t')) newline dataBody];
    tmpFile = [tempname '.tsv'];
    tmpFid = fopen(tmpFile, 'wb');
    fwrite(tmpFid, unicode2native(normalizedBody, 'UTF-8'));
    fclose(tmpFid);
    cleanup = onCleanup(@() delete(tmpFile)); %#ok<NASGU>
    T = readtable(tmpFile, 'FileType', 'text', 'Delimiter', '\t', ...
                  'VariableNamingRule', 'preserve', 'TextType', 'string');
    % UniProt legitimately leaves fields such as EC number, fragment and
    % sequence-modified date empty.  readtable represents those cells as
    % <missing>; normalize them once at the external-data boundary so later
    % filters can safely compare or convert scalar text values.
    T = normalizeTextMissing(T);
end

function T = normalizeTextMissing(T)
    for i = 1:width(T)
        name = T.Properties.VariableNames{i};
        if isstring(T.(name))
            T.(name)(ismissing(T.(name))) = "";
        end
    end
end

function key = computeCacheKey(query, fields)
    payload = [query char(0) fields];
    md = java.security.MessageDigest.getInstance('SHA-256');
    hash = md.digest(uint8(payload));
    key = lower(reshape(dec2hex(hash), 1, []));
end

function rawRows = loadFromCache(path)
    try
        rawRows = readtable(path, 'FileType', 'text', 'Delimiter', '\t', ...
                            'VariableNamingRule', 'preserve', 'TextType', 'string');
        rawRows = normalizeTextMissing(rawRows);
        required = {'accession','organism_name','ec','sequence','length'};
        if ~all(ismember(required, rawRows.Properties.VariableNames)) || ...
                any(strlength(string(rawRows.accession)) > 20)
            rawRows = table();
        end
    catch
        rawRows = table();
    end
end

function saveToCache(path, rawRows)
    parent = fileparts(path);
    if ~exist(parent, 'dir'), mkdir(parent); end
    % A second MATLAB process may query the same key. Write a complete TSV
    % under a unique name and only then replace the shared cache entry, so
    % readers never observe a partially written table.
    tempPath = [tempname(parent) '.tsv'];
    cleanup = onCleanup(@()deleteIfPresent(tempPath)); %#ok<NASGU>
    writetable(rawRows, tempPath, 'FileType', 'text', 'Delimiter', '\t');
    [ok, message] = movefile(tempPath, path, 'f');
    if ~ok, error('fetchCrossSpeciesEnzymes:cacheWrite', '%s', message); end
end

function atomicWriteTable(T, targetPath)
    parent = fileparts(targetPath);
    if ~exist(parent, 'dir'), mkdir(parent); end
    tempPath = [tempname(parent) '.csv'];
    cleanup = onCleanup(@()deleteIfPresent(tempPath)); %#ok<NASGU>
    writetable(T, tempPath);
    [ok, message] = movefile(tempPath, targetPath, 'f');
    if ~ok, error('fetchCrossSpeciesEnzymes:checkpointWrite', '%s', message); end
end

function atomicSaveState(targetPath, candidatesAll, processedKeys, seedKeys)
    parent = fileparts(targetPath);
    if ~exist(parent, 'dir'), mkdir(parent); end
    tempPath = [tempname(parent) '.mat'];
    cleanup = onCleanup(@()deleteIfPresent(tempPath)); %#ok<NASGU>
    save(tempPath, 'candidatesAll', 'processedKeys', 'seedKeys', '-v7');
    [ok, message] = movefile(tempPath, targetPath, 'f');
    if ~ok, error('fetchCrossSpeciesEnzymes:checkpointWrite', '%s', message); end
end
