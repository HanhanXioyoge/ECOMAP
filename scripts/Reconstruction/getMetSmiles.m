function [model, noSMILES] = getMetSmiles(model, parameters)
% getMetSmiles  Retrieve SMILES for metabolites.
% STEP 1: query PubChem by unique metabolite names (with local cache smilesDB.tsv).
% STEP 2: for remaining misses, use MetaNetX IDs (metanetx.chemical) and/or ChEBI IDs
%         to lookup SMILES from local chem_prop.tsv (cols 1=MNX_ID, 3=reference, 9=SMILES).
%
% Inputs
%   model: struct with fields
%       - metNames    : Nx1 cellstr
%       - metMiriams  : Nx1 cell of struct with fields .name and .value
%       - (optional) metSmiles : Nx1 cellstr (will be updated/filled)
%   parameters: struct with field
%       - dataDir     : folder for smilesDB.tsv and chem_prop.tsv
%
% Outputs
%   model   : model.metSmiles updated (keeps existing non-empty entries)
%   noSMILES: cellstr of unique metabolite names that still have no SMILES

    % ---- parameters / paths ----
    if nargin < 2 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.');
        end
    end
    dataDir = parameters.dataDir;
    if ~exist(dataDir, 'dir'), mkdir(dataDir); end
    smilesDBfile  = fullfile(dataDir, 'smilesDB.tsv');
    chemPropPath  = fullfile(findECOMAProot, 'scripts\database\chem_prop.tsv'); % pre-downloaded

    % ---- STEP 1: primary by names (PubChem with local cache) ----
    % Clean names: remove 'prot_' prefix only, then trim spaces
    names = model.metNames;
    if isstring(names), names = cellstr(names); end
    names = regexprep(names, '^prot_.*', '');
    names = cellfun(@strtrim, names, 'UniformOutput', false);

    [uniqueNames, ~, uniqueIdx] = unique(names, 'stable');
    numUnique    = numel(uniqueNames);
    uniqueSmiles = repmat({''}, numUnique, 1);

    % mark empty names as already "matched" (no query)
    isEmptyName = cellfun(@isempty, uniqueNames);
    metMatch = isEmptyName;

    % load local DB (if present)
    if exist(smilesDBfile, 'file') == 2 && dir(smilesDBfile).bytes > 0
        fid = fopen(smilesDBfile, 'r');
        raw = textscan(fid, '%s %s', 'Delimiter', '\t');
        fclose(fid);
        if ~isempty(raw) && numel(raw) >= 2
            db.names  = raw{1};
            db.smiles = raw{2};
            metMatchDB = ismember(uniqueNames, db.names);
            [~, metIdx] = ismember(uniqueNames, db.names);
            hit = metMatchDB & metIdx>0;
            uniqueSmiles(hit) = db.smiles(metIdx(hit));
            metMatch = metMatch | metMatchDB;  % do NOT overwrite the empty-name mask
            fprintf('Local SMILES database loaded (%d entries).\n', numel(db.names));
        else
            fprintf('Local SMILES database found but empty; continuing.\n');
        end
    else
        fprintf('Local SMILES database not found; it will be created/extended at: %s\n', smilesDBfile);
    end

    % query PubChem for remaining names
    opts = weboptions('Timeout', 30, 'ContentType', 'text', ...
                      'CharacterEncoding', 'UTF-8', ...
                      'HeaderFields', {'Accept-Encoding', 'identity'});
    MaxRetries = 3;

    newRecords = {};  % to append once at end: {name, smiles}
    for i = 1:numUnique
        if metMatch(i), continue; end
        q = uniqueNames{i};
        if isempty(q), continue; end

        retry = 0; success = false; sm = '';
        while retry < MaxRetries && ~success
            try
                url  = ['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' q '/property/CanonicalSMILES/TXT'];
                resp = webread(url, opts);
                tok  = regexp(resp, '(^\S+)', 'tokens', 'once');
                sm   = iff(~isempty(tok), tok{1}, '');
                success = true;
            catch ME
                if any(strcmp(ME.identifier, { ...
                        'MATLAB:webservices:HTTP400StatusCodeError', ...
                        'MATLAB:webservices:HTTP404StatusCodeError', ...
                        'MATLAB:webservices:HTTP429StatusCodeError', ...
                        'MATLAB:webservices:HTTP500StatusCodeError'}))
                    success = true; sm = '';
                else
                    retry = retry + 1;
                    pause(0.4 * retry);
                end
            end
        end

        uniqueSmiles{i} = sm;
        if ~isempty(sm)
            newRecords(end+1, :) = {q, sm}; %#ok<AGROW>
            fprintf('Retrieved SMILES from PubChem for %s: %s\n', q, sm);
        end
    end

    % append new PubChem records (if any)
    append_to_db(smilesDBfile, newRecords);

    % ---- STEP 2: secondary via MetaNetX chem_prop.tsv (MNX_ID / ChEBI → SMILES) ----
    if any(cellfun(@isempty, uniqueSmiles))
        if ~isfile(chemPropPath)
            warning('chem_prop.tsv not found at %s. Skipping MetaNetX step.', chemPropPath);
        else
            % Build two maps: MNX_ID->SMILES and chebi:NNNNN->SMILES (lowercase)
            [mnxMap, chebiMap] = build_mnx_and_chebi_maps(chemPropPath);

            newRecords2 = {};
            for i = 1:numUnique
                if ~isempty(uniqueSmiles{i}), continue; end

                % first index of this unique name in original list
                idxList = find(strcmp(names, uniqueNames{i}));
                if isempty(idxList), continue; end
                firstIdx = idxList(1);

                if firstIdx <= numel(model.metMiriams) && ~isempty(model.metMiriams{firstIdx})
                    miriam = model.metMiriams{firstIdx};

                    % gather candidate MNX and ChEBI IDs (multiple allowed, exact match only)
                    candMnx  = extract_mnx_ids_from_miriam(miriam);             % e.g., {'MNXM4939',...}
                    candChEB = extract_chebi_ids_from_miriam(miriam, true);     % e.g., {'chebi:62958',...} (lowercase)

                    % 1) try MNX_IDs in order
                    sm = '';
                    for k = 1:numel(candMnx)
                        id = candMnx{k};
                        if isKey(mnxMap, id)
                            sm = mnxMap(id);
                            if ~isempty(sm), break; end
                        end
                    end

                    % 2) if still empty, try ChEBI IDs in order
                    if isempty(sm)
                        for k = 1:numel(candChEB)
                            key = candChEB{k};     % already 'chebi:NNNNN' lowercase
                            if isKey(chebiMap, key)
                                sm = chebiMap(key);
                                if ~isempty(sm), break; end
                            end
                        end
                    end

                    % fill if found
                    if ~isempty(sm)
                        uniqueSmiles{i} = sm;
                        if ~isempty(uniqueNames{i})
                            newRecords2(end+1,:) = {uniqueNames{i}, sm}; %#ok<AGROW>
                        end
                        fprintf('MetaNetX: Retrieved SMILES for %s.\n', uniqueNames{i});
                    end
                end
            end

            % append new MetaNetX records
            append_to_db(smilesDBfile, newRecords2);
        end
    end

%% STEP 3: MetaNetX chem_info JSON-LD fallback (fetch web page and parse "smiles")
    % Only try for the remaining unresolved unique names.
    if any(cellfun(@isempty, uniqueSmiles))
        newRecords3 = {};
        MaxWebRetries = 3;
    
        for i = 1:numUnique
            if ~isempty(uniqueSmiles{i}), continue; end
    
            % find first occurrence of this unique name in original list
            idxList = find(strcmp(names, uniqueNames{i}));
            if isempty(idxList), continue; end
            firstIdx = idxList(1);
    
            if firstIdx <= numel(model.metMiriams) && ~isempty(model.metMiriams{firstIdx})
                miriam = model.metMiriams{firstIdx};
    
                % gather candidate MNX IDs (multiple allowed, exact pattern 'MNXM\d+')
                candMnx = extract_mnx_ids_from_miriam(miriam);
                if isempty(candMnx), continue; end
    
                % try each MNX in order until one returns a SMILES
                sm = '';
                for k = 1:numel(candMnx)
                    mnxID = candMnx{k};
                    % polite retry loop per MNX page
                    for r = 1:MaxWebRetries
                        try
                            sm = fetch_smiles_from_mnx_page(mnxID);
                            break; % break retry loop
                        catch
                            sm = '';
                            pause(0.25*r); % mild backoff
                        end
                    end
                    if ~isempty(sm), break; end % break MNX candidate loop
                end
    
                if ~isempty(sm)
                    uniqueSmiles{i} = sm;
                    if ~isempty(uniqueNames{i})
                        newRecords3(end+1,:) = {uniqueNames{i}, sm}; %#ok<AGROW>
                    end
                    fprintf('MetaNetX JSON-LD: Retrieved SMILES for %s.\n', uniqueNames{i});
                end
            end
        end
    
        % append newly found records in one I/O
        append_to_db(smilesDBfile, newRecords3);
    end

    % ---- reassemble outputs ----
    newSmiles = uniqueSmiles(uniqueIdx);
    noSMILES  = uniqueNames(cellfun(@isempty, uniqueSmiles));

    % fill model.metSmiles (keep existing non-empty)
    if ~isfield(model, 'metSmiles') || isempty(model.metSmiles)
        model.metSmiles = newSmiles;
    else
        ms = model.metSmiles;
        if isstring(ms), ms = cellstr(ms); end
        ms = ms(:);
        if numel(ms) ~= numel(newSmiles), ms = repmat({''}, numel(newSmiles), 1); end
        emptyMask = cellfun(@isempty, ms);
        ms(emptyMask) = newSmiles(emptyMask);
        model.metSmiles = ms;
    end

    % simple report (exclude empty-name rows from denominator)
    denom = sum(~isEmptyName);
    successRatio = iff(denom==0, 0, sum(~cellfun(@isempty, uniqueSmiles(~isEmptyName))) / denom);
    fprintf('SMILES found for %.1f%% of unique non-empty metabolite names.\n', 100*successRatio);
end

% =================== local helpers ===================
function v = iff(cond, a, b)
    if cond, v = a; else, v = b; end
end

function append_to_db(dbfile, records)
% Append {name,smiles} rows to smilesDB.tsv in one shot; ignore empty input.
    if isempty(records), return; end
    fid = fopen(dbfile, 'a');
    c = onCleanup(@() fclose(fid));
    for k = 1:size(records,1)
        if ~isempty(records{k,1}) && ~isempty(records{k,2})
            fprintf(fid, '%s\t%s\n', records{k,1}, records{k,2});
        end
    end
end

function arr = field_to_cellstr(s, fname)
% Convert s.(fname) into cellstr, handling char/string/cell.
    arr = {};
    if ~isfield(s, fname), return; end
    v = s.(fname);
    if iscell(v)
        for k = 1:numel(v), arr{k,1} = char(v{k}); end
    else
        arr = {char(v)};
    end
end

function ids = extract_mnx_ids_from_miriam(miriam)
% Collect *all* MNX IDs from miriam where name suggests MetaNetX entry first,
% then fallback to scan all values. Exact pattern 'MNXM\d+'.
    ids = {};
    mNames = field_to_cellstr(miriam, 'name');
    mVals  = field_to_cellstr(miriam, 'value');

    % priority 1: values paired with names that contain 'metanetx.chemical'
    pri = find(contains(lower(mNames), 'metanetx.chemical'));
    for p = pri(:).'
        tok = regexp(mVals{p}, 'MNXM\d+', 'match');
        if ~isempty(tok), ids = [ids, tok]; end%#ok<AGROW>
    end

    % fallback: scan all values
    if isempty(ids)
        for j = 1:numel(mVals)
            tok = regexp(mVals{j}, 'MNXM\d+', 'match');
            if ~isempty(tok), ids = [ids, tok]; end %#ok<AGROW>
        end
    end

    if ~isempty(ids)
        ids = unique(ids, 'stable');
    end
end

function ids = extract_chebi_ids_from_miriam(miriam, toLower)
% Collect *all* ChEBI IDs from miriam (exact match), accept variants:
% 'CHEBI:12345', 'chebi:12345', or raw number in a value paired with name containing 'chebi'.
% Return as 'chebi:12345' if toLower==true, otherwise 'CHEBI:12345'.
    ids = {};
    mNames = field_to_cellstr(miriam, 'name');
    mVals  = field_to_cellstr(miriam, 'value');

    norm = @(s) regexprep(char(s), '^\s*(?i)chebi:\s*', '', 'once'); % remove CHEBI: prefix (any case)
    emit = @(num) iff(toLower, ['chebi:' num], ['CHEBI:' num]);

    % priority 1: values paired with names containing 'chebi'
    pri = find(contains(lower(mNames), 'chebi'));
    for p = pri(:).'
        s = mVals{p};
        % 1) direct CHEBI:NNNNN
        tok = regexp(s, '(?i)CHEBI:\s*(\d+)', 'tokens');
        for t = 1:numel(tok), ids{end+1} = emit(tok{t}{1}); end %#ok<AGROW>
        % 2) bare number if the name says 'chebi'
        tok2 = regexp(s, '^\s*(\d+)\s*$', 'tokens', 'once');
        if ~isempty(tok2), ids{end+1} = emit(tok2{1}); end %#ok<AGROW>
    end

    % fallback: scan all values for explicit CHEBI:NNNNN
    for j = 1:numel(mVals)
        tok = regexp(mVals{j}, '(?i)CHEBI:\s*(\d+)', 'tokens');
        for t = 1:numel(tok), ids{end+1} = emit(tok{t}{1}); end %#ok<AGROW>
    end

    if ~isempty(ids)
        ids = unique(ids, 'stable');
    end
end

function [mnxMap, chebiMap] = build_mnx_and_chebi_maps(chemPropPath)
% Build two containers.Map from chem_prop.tsv:
%   - mnxMap  : MNX_ID (e.g., 'MNXM4939') -> SMILES
%   - chebiMap: 'chebi:NNNNN' (lowercase) -> SMILES
% Stream-read only Var1 (MNX_ID), Var3 (reference), Var9 (SMILES).

    % allow .tsv extension; TreatAsMissing contains only non-empty tokens
    try
        ds = tabularTextDatastore(chemPropPath, ...
            'Delimiter', '\t', ...
            'TextType', 'string', ...
            'ReadVariableNames', false, ...
            'NumHeaderLines', 0, ...
            'TreatAsMissing', {'NA'}, ...
            'CommentStyle', '#', ...
            'FileExtensions', {'.tsv','.txt','.tab'});
    catch
        % older MATLAB compatibility
        ds = datastore(chemPropPath, 'Type', 'tabulartext', ...
            'Delimiter', '\t', ...
            'TextType', 'string', ...
            'ReadVariableNames', false, ...
            'NumHeaderLines', 0, ...
            'TreatAsMissing', {'NA'}, ...
            'CommentStyle', '#', ...
            'FileExtensions', {'.tsv','.txt','.tab'});
    end

    vars = ds.VariableNames;  % {'Var1','Var2',...}
    if numel(vars) < 9
        error('chem_prop.tsv appears to have fewer than 9 columns.');
    end
    ds.SelectedVariableNames = vars([1 3 9]); % MNX_ID, reference, SMILES

    mnxMap   = containers.Map('KeyType','char','ValueType','char');
    chebiMap = containers.Map('KeyType','char','ValueType','char');

    while hasdata(ds)
        T = read(ds);  % T.(vars{1})=MNX_ID, T.(vars{3})=reference, T.(vars{9})=SMILES
        if isempty(T), continue; end

        idCol = T.(vars{1});
        refCol= T.(vars{3});
        smCol = T.(vars{9});

        % keep only rows with non-empty MNX_ID and SMILES
        ok = (strlength(idCol) > 0) & (strlength(smCol) > 0);
        if ~any(ok), continue; end

        idCol = idCol(ok);
        refCol= refCol(ok);
        smCol = smCol(ok);

        % fill MNX_ID -> SMILES
        for k = 1:numel(idCol)
            key = char(idCol(k));
            if ~isKey(mnxMap, key)
                mnxMap(key) = char(smCol(k));
            end
        end

        % from reference column, extract chebi:NNNNN (lowercase), map to SMILES
        for k = 1:numel(refCol)
            r = lower(char(refCol(k)));
            if isempty(r), continue; end
            toks = regexp(r, 'chebi:\d+', 'match');  % may return multiple tokens per row
            if isempty(toks), continue; end
            sm = char(smCol(k));
            for t = 1:numel(toks)
                ck = toks{t}; % already lowercase
                if ~isKey(chebiMap, ck)
                    chebiMap(ck) = sm;
                end
            end
        end
    end
end

function sm = fetch_smiles_from_mnx_page(mnxID)
% Fetch "smiles" from the JSON-LD block of https://www.metanetx.org/chem_info/<MNXID>
% Returns '' if not found.

    % normalize ID like 'MNXM146479'
    mnxID = upper(char(mnxID));
    if isempty(regexp(mnxID,'^MNXM\d+$','once'))
        error('fetch_smiles_from_mnx_page:BadID', 'Not a valid MNX ID: %s', mnxID);
    end

    url = ['https://www.metanetx.org/chem_info/' mnxID];
    opts = weboptions('Timeout', 30, ...
                      'ContentType', 'text', ...
                      'CharacterEncoding', 'UTF-8', ...
                      'HeaderFields', { ...
                         'Accept','text/html,application/xhtml+xml,application/json;q=0.9,*/*;q=0.8'; ...
                         'Accept-Encoding','identity'; ...
                         'User-Agent','Mozilla/5.0 (MATLAB)'} );

    html = webread(url, opts);

    % Extract the first JSON-LD <script type="application/ld+json"> ... </script>
    % Use (?s) to let '.' match newlines
    toks = regexp(html, '(?s)<script[^>]*type="application/ld\+json"[^>]*>\s*(\{.*?\})\s*</script>', 'tokens');

    sm = '';
    if isempty(toks), return; end

    % Prefer a block that has "smiles" or @type MolecularEntity
    for t = 1:numel(toks)
        txt = strtrim(toks{t}{1});
        % Some pages may contain HTML entities; jsondecode usually tolerates them as plain text
        try
            j = jsondecode(txt);
        catch
            continue
        end
        if isstruct(j)
            if isfield(j,'smiles')
                sm = coerce_smiles_field(j.smiles);
                if ~isempty(sm), return; end
            end
            % some pages could wrap in arrays; handle minimal cases
        elseif iscell(j)
            for c = 1:numel(j)
                if isstruct(j{c}) && isfield(j{c},'smiles')
                    sm = coerce_smiles_field(j{c}.smiles);
                    if ~isempty(sm), return; end
                end
            end
        end
    end
end

function sm = coerce_smiles_field(val)
% Normalize "smiles" JSON value to a single char row (take first if array)
    if ischar(val) || isstring(val)
        sm = char(val);
    elseif iscell(val) && ~isempty(val)
        sm = char(val{1});
    else
        sm = '';
    end
    sm = strtrim(sm);
end
