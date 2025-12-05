function [model, noSMILES, noInChIKey] = getMetinfo(model, verbosity, methods, parameters)
% getMetinfo
% Populate per-metabolite SMILES and InChIKey using up to three sources:
%   STEP A) MetaNetX chem_prop.tsv (via MNX_ID / ChEBI cross-refs)
%   STEP B) ChEBI public API by CHEBI:ID parsed from model.metMiriams
%   STEP C) PubChem PUG-REST TXT endpoints by metabolite *names*
%
% SYNTAX
%   [model, noSMILES, noInChIKey] = getMetinfo(model, verbosity, parameters)
%   [model, noSMILES, noInChIKey] = getMetinfo(model, verbosity, parameters, methods)
%
% INPUTS
%   model       : struct with at least
%                 - metNames    : Nx1 cellstr/strings
%                 - metMiriams  : Nx1 cell of struct(s) with .name / .value
%                 - (optional) metSmiles   : Nx1 cellstr (kept if non-empty)
%                 - (optional) metInChIKey : Nx1 cellstr (kept if non-empty)
%   verbosity   : 0=silent, 1=info, 2=debug (controls prints; NOT part of 'parameters').
%   parameters  : struct with at least
%                 - dataDir : folder path for the local cache 'metInfo.tsv'
%   methods     : which steps to run. Accepts:
%                 - string like 'ABC', 'AC', 'B' (case-insensitive). Default 'ABC'.
%                 - logical/numeric vector [A B C], e.g. [1 0 1].
%                 - struct with fields .A .B .C (logical).
%
% OUTPUTS
%   model       : with model.metSmiles and model.metInChIKey filled (existing non-empties preserved)
%   noSMILES    : unique names still missing SMILES after selected steps
%   noInChIKey  : unique names still missing InChIKey after selected steps
%
% NOTES
%   - Local cache file: <dataDir>/metInfo.tsv  (columns: name \t smiles \t inchikey).
%     It is read once and re-written atomically at the end if updated.
%   - Name cleaning: drop names that match '^prot_.*' entirely, then strtrim.
%   - STEP B endpoint:
%       https://www.ebi.ac.uk/chebi/backend/api/public/compound/CHEBI:12345/
%   - STEP C endpoints (TXT):
%       /compound/name/<NAME-encoded>/property/CanonicalSMILES/TXT
%       /compound/name/<NAME-encoded>/property/InChIKey/TXT
%     We send ?name_type=word to improve recall; path segment is %20-encoded (NOT '+').
%   - All InChIKeys are normalized by stripping a leading 'InChIKey=' (case-insensitive).

    % -------- defaults / guards --------
    if nargin < 2 || isempty(verbosity),  verbosity = 1; end
    
    if (nargin == 3 || (nargin >= 3 && isempty(methods))) && isstruct(methods) ...
            && ~isfield(methods,'A') && ~isfield(methods,'B') && ~isfield(methods,'C')
        parameters = methods;   % reassign 3rd arg to 'parameters'
        methods    = [];        % mark 'methods' as omitted
    end
    
    if nargin < 4 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end
    
    % Default when 'methods' is omitted: use ONLY local metInfo.tsv cache.
    % That means disabling A (chem_prop.tsv), B (ChEBI API), and C (PubChem).
    if nargin < 3 || isempty(methods), methods = [0 0 0]; end
    [useA, useB, useC] = normalize_methods(methods);

    if ~isfield(parameters, 'dataDir') || ~isfolder(parameters.dataDir)
        if ~isfield(parameters,'dataDir'), error('parameters.dataDir is required.'); end
        mkdir(parameters.dataDir);
    end

    dataDir      = parameters.dataDir;
    metInfoFile  = fullfile(dataDir, 'metInfo.tsv');  % renamed from smilesDB.tsv
    chemPropPath = fullfile(findECOMAProot, 'scripts','database','chem_prop.tsv');

    % -------- names preprocessing --------
    names = model.metNames;
    if isstring(names), names = cellstr(names); end
    names = regexprep(names, '^prot_.*', '');          % remove whole-string that starts with 'prot_'
    names = cellfun(@strtrim, names, 'UniformOutput', false);

    [uniqueNames, ~, uniqueIdx] = unique(names, 'stable');
    numUnique     = numel(uniqueNames);
    isEmptyName   = cellfun(@isempty, uniqueNames);

    uniqueSmiles  = repmat({''}, numUnique, 1);
    uniqueInChI   = repmat({''}, numUnique, 1);
    hasInfo       = isEmptyName;  % mask of names already resolved (empty name = resolved/ignored)

    % -------- load local cache: metInfo.tsv (name,smiles,inchikey) --------
    dbMap = containers.Map('KeyType','char','ValueType','any'); % value: struct('smiles',..,'inchikey',..)
    if isfile(metInfoFile) && dir(metInfoFile).bytes > 0
        fid = fopen(metInfoFile, 'r');
        c = onCleanup(@() fclose(fid));
        raw = textscan(fid, '%s%s%s', 'Delimiter', '\t', 'ReturnOnError', false);
        if isempty(raw{3}), raw{3} = repmat({''}, numel(raw{1}), 1); end
        for i = 1:numel(raw{1})
            nm = orEmpty(raw{1}{i});
            sm = orEmpty(raw{2}{i});
            ik = clean_inchikey(orEmpty(raw{3}{i}));  % normalize prefix here
            if ~isempty(nm)
                dbMap(nm) = struct('smiles', sm, 'inchikey', ik);
            end
        end
        logv(verbosity,1,'Loaded local metInfo DB (%d entries).\n', numel(raw{1}));
    else
        logv(verbosity,1,'Local metInfo DB not found; it will be created at: %s\n', metInfoFile);
    end
    dirty = false; % whether we need to rewrite metInfo.tsv at the end

    % Seed from local DB by name
    for i = 1:numUnique
        if hasInfo(i), continue; end
        nm = uniqueNames{i};
        if isempty(nm), continue; end
        if isKey(dbMap, nm)
            rec = dbMap(nm);
            uniqueSmiles{i} = orEmpty(rec.smiles);
            uniqueInChI{i}  = clean_inchikey(orEmpty(rec.inchikey)); % ensure cleaned
            if verbosity >= 1
                logv(verbosity,1,'[cache] %s -> SMILES:%s InChIKey:%s\n', nm, ...
                    safeShort(uniqueSmiles{i}), safeShort(uniqueInChI{i}));
            end
            hasInfo(i) = (~isempty(uniqueSmiles{i}) && ~isempty(uniqueInChI{i}));
        end
    end

    % =====================================================================
    % STEP A) chem_prop.tsv (MNX_ID / ChEBI → SMILES, InChIKey)
    % =====================================================================
    if useA
        if any(~hasInfo & ~isEmptyName)
            if ~isfile(chemPropPath)
                logv(verbosity,1,'chem_prop.tsv not found at %s. Skipping MetaNetX step.\n', chemPropPath);
            else
                [mnxMap, chebiMap] = build_mnx_chebi_maps(chemPropPath); % maps to struct(smiles,inchikey)

                for i = 1:numUnique
                    if hasInfo(i) || isEmptyName(i), continue; end

                    idxList = find(strcmp(names, uniqueNames{i}));
                    if isempty(idxList), continue; end
                    firstIdx = idxList(1);

                    sm = uniqueSmiles{i};
                    ik = uniqueInChI{i};

                    if firstIdx <= numel(model.metMiriams) && ~isempty(model.metMiriams{firstIdx})
                        miriam = model.metMiriams{firstIdx};

                        candMnx  = extract_mnx_ids_from_miriam(miriam);             % {'MNXM123',...}
                        candChEB = extract_chebi_ids_from_miriam(miriam, true);     % {'chebi:12345',...} lowercase

                        % Try MNX IDs
                        if isempty(sm) || isempty(ik)
                            for k = 1:numel(candMnx)
                                id = candMnx{k};
                                if isKey(mnxMap, id)
                                    rec = mnxMap(id);
                                    if isempty(sm), sm = orEmpty(rec.smiles); end
                                    if isempty(ik), ik = clean_inchikey(orEmpty(rec.inchikey)); end
                                    if ~isempty(sm) && ~isempty(ik), break; end
                                end
                            end
                        end

                        % Try ChEBI IDs in chem_prop
                        if isempty(sm) || isempty(ik)
                            for k = 1:numel(candChEB)
                                key = candChEB{k};
                                if isKey(chebiMap, key)
                                    rec = chebiMap(key);
                                    if isempty(sm), sm = orEmpty(rec.smiles); end
                                    if isempty(ik), ik = clean_inchikey(orEmpty(rec.inchikey)); end
                                    if ~isempty(sm) && ~isempty(ik), break; end
                                end
                            end
                        end
                    end

                    if ~isempty(sm) || ~isempty(ik)
                        uniqueSmiles{i} = sm;
                        uniqueInChI{i}  = ik;
                        if ~isempty(sm) && ~isempty(ik)
                            hasInfo(i) = true;
                        end
                        nm = uniqueNames{i};
                        changed = upsert_metinfo(dbMap, nm, sm, ik, verbosity);
                        dirty = dirty || changed;

                        if verbosity >= 1
                            logv(verbosity,1,'[chem_prop] %s -> SMILES:%s InChIKey:%s\n', nm, safeShort(sm), safeShort(ik));
                        end
                    end
                end
            end
        end
    else
        if verbosity>=1, logv(verbosity,1,'[skip] STEP A disabled by methods.\n'); end
    end

    % =====================================================================
    % STEP B) ChEBI by CHEBI:ID (from MIRIAM), for residues still missing either field
    % =====================================================================
    if useB
        needChEBI = find(~hasInfo & ~isEmptyName);
        if ~isempty(needChEBI)
            CE_BASE = "https://www.ebi.ac.uk/chebi/backend/api/public/compound/";
            CE_UA   = 'ECOMAP-getMetinfo/1.0 (ChEBI resolver)';
            ceOpts  = weboptions('Timeout', 15, ...
                                 'HeaderFields', {'Accept','application/json'; 'User-Agent', CE_UA}, ...
                                 'ContentType','text');  % raw text; decode manually
            MAX_RETRY_PER_ID = 1;
            MIN_INTERVAL_SEC = 0.15;   % be polite
            lastTickCE = [];

            for idx = needChEBI(:).'
                nm = uniqueNames{idx};
                if isempty(nm), continue; end

                idxList = find(strcmp(names, nm));
                firstIdx = idxList(1);
                if firstIdx > numel(model.metMiriams) || isempty(model.metMiriams{firstIdx})
                    continue;
                end
                miriam = model.metMiriams{firstIdx};

                candChEBI = extract_chebi_ids_from_miriam(miriam, false);  % {'CHEBI:12345',...}
                if isempty(candChEBI), continue; end

                sm = uniqueSmiles{idx};
                ik = uniqueInChI{idx};
                if ~isempty(sm) && ~isempty(ik)
                    continue;
                end

                got = false;
                for k = 1:numel(candChEBI)
                    chebiId = candChEBI{k};
                    if isempty(chebiId), continue; end

                    lastTickCE = rate_limit(lastTickCE, MIN_INTERVAL_SEC);
                    url = CE_BASE + encode_chebi_for_path(chebiId) + "/";

                    hit = false; lastErr = "";
                    for r = 0:MAX_RETRY_PER_ID
                        try
                            txt = webread(url, ceOpts);          % raw text
                            data = safe_jsondecode(txt);         % may be []
                            [s2, i2] = extract_smiles_inchikey(data, txt);
                            if isempty(sm) && ~isempty(s2), sm = s2; end
                            if isempty(ik) && ~isempty(i2), ik = clean_inchikey(i2); end
                            hit = (~isempty(sm) || ~isempty(ik));
                            break;
                        catch ME
                            lastErr = ME.message;
                            pause(0.4*(r+1));
                        end
                    end

                    if hit
                        got = true;
                        if verbosity >= 1
                            logv(verbosity,1,'[chebi] %s (%s) -> SMILES:%s InChIKey:%s\n', nm, chebiId, safeShort(sm), safeShort(ik));
                        end
                        break;
                    else
                        if verbosity >= 2
                            logv(verbosity,2,'[chebi][miss] %s (%s) : %s\n', nm, chebiId, lastErr);
                        end
                    end
                end

                if got
                    uniqueSmiles{idx} = sm;
                    uniqueInChI{idx}  = ik;
                    hasInfo(idx) = (~isempty(sm) && ~isempty(ik));

                    changed = upsert_metinfo(dbMap, nm, sm, ik, verbosity);
                    dirty = dirty || changed;
                end
            end
        end
    else
        if verbosity>=1, logv(verbosity,1,'[skip] STEP B disabled by methods.\n'); end
    end

    % =====================================================================
    % STEP C) PubChem by metabolite names (TXT endpoints, split properties)
    %         - CanonicalSMILES:  /property/CanonicalSMILES/TXT
    %         - InChIKey       :  /property/InChIKey/TXT
    %         - name_type=word to improve recall
    %         - 503/429: exponential backoff; early abort on too many other failures
    % =====================================================================
    if useC
        needPubChem = find(~hasInfo & ~isEmptyName);
        if ~isempty(needPubChem)
            PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/";
            pcTxtOpts = weboptions('Timeout', 30);
            MAX_RETRY_PER_MET = 10;       % retry count for non-404 errors
            MIN_INTERVAL_SEC  = 2;     % polite pacing
            MAX_CONSEC_FAILS  = 10;      % early abort threshold

            lastTickPC  = [];
            failStreak  = 0;
            abortPC     = false;

            logv(verbosity,1,'[pubchem/txt] pending = %d\n', numel(needPubChem));

            for idx = needPubChem(:).'
                if abortPC, break; end

                nm = uniqueNames{idx};
                if isempty(nm), continue; end

                sm = uniqueSmiles{idx};
                ik = uniqueInChI{idx};

                % CanonicalSMILES (only if still missing)
                if isempty(sm)
                    lastTickPC = rate_limit(lastTickPC, MIN_INTERVAL_SEC);
                    qname = nm;
                    urlS  = PUG_BASE + qname + "/property/CanonicalSMILES/TXT";

                    okS = false; lastErrS = "";
                    for r = 0:MAX_RETRY_PER_MET
                        try
                            txtS = webread(urlS, pcTxtOpts);
                            % sm2  = parse_pubchem_txt_property(txtS, 'CanonicalSMILES');
                            sm2 = regexp(txtS,'(^\S*)\n','once','tokens');
                            if ~isempty(sm2), sm = sm2{1,1}; end
                            failStreak = 0;
                            okS = ~isempty(sm);
                            break;
                        catch ME
                            lastErrS = ME.message;
                            if contains(lastErrS,'404') || contains(lower(lastErrS),'not found')
                                okS = false; break;
                            end
                            if contains(lastErrS,'503') || contains(lastErrS,'429')
                                pause(min(5, 0.5*(2^r)));
                            else
                                failStreak = failStreak + 1;
                                if failStreak >= MAX_CONSEC_FAILS
                                    logv(verbosity,1,'[pubchem/txt] too many consecutive failures (%d). Aborting PubChem step early.\n', failStreak);
                                    abortPC = true; break;
                                end
                                pause(0.3*(r+1));
                            end
                        end
                    end
                    if verbosity >= 2 && ~okS
                        logv(verbosity,2,'[pubchem/txt][miss] CanonicalSMILES %s : %s\n', nm, lastErrS);
                    end
                end

                % InChIKey (only if still missing)
                if ~abortPC && isempty(ik)
                    lastTickPC = rate_limit(lastTickPC, MIN_INTERVAL_SEC);
                    qname = nm;
                    urlK  = PUG_BASE + qname + "/property/InChIKey/TXT";

                    okK = false; lastErrK = "";
                    for r = 0:MAX_RETRY_PER_MET
                        try
                            txtK = webread(urlK, pcTxtOpts);
                            % ik2  = parse_pubchem_txt_property(txtK, 'InChIKey');
                            ik2 = regexp(txtK,'(^\S*)\n','once','tokens');
                            if ~isempty(ik2), ik = ik2{1,1}; end
                            failStreak = 0;
                            okK = ~isempty(ik);
                            break;
                        catch ME
                            lastErrK = ME.message;
                            if contains(lastErrK,'404') || contains(lower(lastErrK),'not found')
                                okK = false; break;
                            end
                            if contains(lastErrK,'503') || contains(lastErrK,'429')
                                pause(min(5, 0.5*(2^r)));
                            else
                                failStreak = failStreak + 1;
                                if failStreak >= MAX_CONSEC_FAILS
                                    logv(verbosity,1,'[pubchem/txt] too many consecutive failures (%d). Aborting PubChem step early.\n', failStreak);
                                    abortPC = true; break;
                                end
                                pause(0.3*(r+1));
                            end
                        end
                    end
                    if verbosity >= 2 && ~okK
                        logv(verbosity,2,'[pubchem/txt][miss] InChIKey %s : %s\n', nm, lastErrK);
                    end
                end

                % Write back if anything was obtained
                if ~abortPC && (~isempty(sm) || ~isempty(ik))
                    uniqueSmiles{idx} = sm;
                    uniqueInChI{idx}  = ik;
                    if ~isempty(sm) && ~isempty(ik), hasInfo(idx) = true; end

                    changed = upsert_metinfo(dbMap, nm, sm, ik, verbosity);
                    dirty = dirty || changed;

                    if verbosity >= 1
                        logv(verbosity,1,'[pubchem/txt] %s -> SMILES:%s InChIKey:%s\n', nm, safeShort(sm), safeShort(ik));
                    end
                end
            end
        end
    else
        if verbosity>=1, logv(verbosity,1,'[skip] STEP C disabled by methods.\n'); end
    end

    % -------- persist local cache once (avoid duplicates) --------
    if dirty
        save_metinfo(dbMap, metInfoFile);
        logv(verbosity,1,'Saved metInfo DB (%d names) to %s\n', dbMap.Count, metInfoFile);
    end

    % -------- assemble outputs back to model fields --------
    newSmiles = uniqueSmiles(uniqueIdx);
    newInChI  = uniqueInChI(uniqueIdx);

    % metSmiles
    if ~isfield(model, 'metSmiles') || isempty(model.metSmiles)
        model.metSmiles = newSmiles;
    else
        ms = model.metSmiles; if isstring(ms), ms = cellstr(ms); end; ms = ms(:);
        if numel(ms) ~= numel(newSmiles), ms = repmat({''}, numel(newSmiles), 1); end
        mask = cellfun(@isempty, ms); ms(mask) = newSmiles(mask);
        model.metSmiles = ms;
    end

    % metInChIKey
    if ~isfield(model, 'metInChIKey') || isempty(model.metInChIKey)
        model.metInChIKey = cellfun(@clean_inchikey, newInChI, 'UniformOutput', false);
    else
        mk = model.metInChIKey; if isstring(mk), mk = cellstr(mk); end; mk = mk(:);
        if numel(mk) ~= numel(newInChI), mk = repmat({''}, numel(newInChI), 1); end
        mask = cellfun(@isempty, mk);
        cleaned = cellfun(@clean_inchikey, newInChI(mask), 'UniformOutput', false);
        mk(mask) = cleaned;
        mk(~mask) = cellfun(@clean_inchikey, mk(~mask), 'UniformOutput', false);
        model.metInChIKey = mk;
    end

    % -------- coverage (include partial hits) --------
    denom = sum(~isEmptyName);
    haveS  = ~cellfun(@isempty, uniqueSmiles(~isEmptyName));
    haveI  = ~cellfun(@isempty, uniqueInChI(~isEmptyName));
    haveAny  = haveS | haveI;
    haveBoth = haveS & haveI;

    covS    = ifthen(denom==0, 0, sum(haveS)    / denom);
    covI    = ifthen(denom==0, 0, sum(haveI)    / denom);
    covAny  = ifthen(denom==0, 0, sum(haveAny)  / denom);
    covBoth = ifthen(denom==0, 0, sum(haveBoth) / denom);

    if verbosity >= 1
        logv(verbosity,1,'Coverage (unique non-empty names):\n');
        logv(verbosity,1,'  SMILES:    %.1f%%\n', 100*covS);
        logv(verbosity,1,'  InChIKey:  %.1f%%\n', 100*covI);
        logv(verbosity,1,'  Any(one):  %.1f%%\n', 100*covAny);
        logv(verbosity,1,'  Both:      %.1f%%\n', 100*covBoth);
    end

    % report lists
    noSMILES   = uniqueNames(~isEmptyName & cellfun(@isempty, uniqueSmiles));
    noInChIKey = uniqueNames(~isEmptyName & cellfun(@isempty, uniqueInChI));
end

% =================== helpers ===================

function [A,B,C] = normalize_methods(methods)
% Normalize 'methods' selector into three booleans A/B/C.
    A = true; B = true; C = true; % default 'ABC'
    if ischar(methods) || isstring(methods)
        s = upper(char(methods));
        A = contains(s,'A'); B = contains(s,'B'); C = contains(s,'C');
    elseif isnumeric(methods) || islogical(methods)
        v = logical(methods(:).'); v(end+1:3) = false;  % pad
        A = v(1); B = v(2); C = v(3);
    elseif isstruct(methods)
        if isfield(methods,'A'), A = logical(methods.A); end
        if isfield(methods,'B'), B = logical(methods.B); end
        if isfield(methods,'C'), C = logical(methods.C); end
    end
end

function enc = encode_path_segment(str)
% Encode a string for a URL *path segment*:
% - Spaces become %20 (NOT '+'); '+' in the name becomes %2B.
% - Leave %XX sequences as-is.
    enc = char(java.net.URLEncoder.encode(char(string(str)), 'UTF-8'));
    enc = strrep(enc, '+', '%20');
end

function val = parse_pubchem_txt_property(txt, propName)
% Parse PUG-REST /TXT outputs for property tables.
% Handles:
%   - header + tab-separated rows (CID \t <prop>)
%   - single-line outputs (just value or "Prop: value")
%   - regex fallback for InChIKey pattern
    val = '';
    if isempty(txt), return; end
    lines = regexp(txt, '\r?\n', 'split');
    lines = lines(~cellfun(@isempty, strtrim(lines)));
    if isempty(lines), return; end

    hdr = regexp(lines{1}, '\t', 'split');
    if numel(hdr) > 1 && any(strcmpi(strtrim(hdr), propName))
        colNames = strtrim(hdr);
        j = find(strcmpi(colNames, propName), 1, 'first');
        if ~isempty(j)
            for i = 2:numel(lines)
                toks = regexp(lines{i}, '\t', 'split');
                if numel(toks) >= j
                    cand = strtrim(toks{j});
                    if ~isempty(cand), val = cand; return; end
                end
            end
        end
    end

    if numel(lines) == 1
        m = regexp(lines{1}, ['^\s*' propName '\s*:\s*(.+)$'], 'tokens', 'once');
        if ~isempty(m), val = strtrim(m{1}); if ~isempty(val), return; end
        else, val = strtrim(lines{1}); if ~isempty(val), return; end
        end
    end

    if strcmpi(propName,'InChIKey')
        m = regexp(strjoin(lines,' '), '([A-Z]{14}-[A-Z]{10}-[A-Z])', 'tokens', 'once');
        if ~isempty(m), val = m{1}; return; end
    end
end

function [mnxMap, chebiMap] = build_mnx_chebi_maps(chemPropPath)
% Stream-read chem_prop.tsv to build:
%   mnxMap   : 'MNXM123' -> struct('smiles',..., 'inchikey',...)
%   chebiMap : 'chebi:12345' (lowercase) -> struct('smiles',..., 'inchikey',...)
% Uses columns: 1=MNX_ID, 3=reference, 8=InChIKey, 9=SMILES
    try
        ds = tabularTextDatastore(chemPropPath, ...
            'Delimiter', '\t', 'TextType', 'string', ...
            'ReadVariableNames', false, 'NumHeaderLines', 0, ...
            'TreatAsMissing', {'NA'}, 'CommentStyle', '#', ...
            'FileExtensions', {'.tsv','.txt','.tab'});
    catch
        ds = datastore(chemPropPath, 'Type', 'tabulartext', ...
            'Delimiter', '\t', 'TextType', 'string', ...
            'ReadVariableNames', false, 'NumHeaderLines', 0, ...
            'TreatAsMissing', {'NA'}, 'CommentStyle', '#', ...
            'FileExtensions', {'.tsv','.txt','.tab'});
    end
    allVars = ds.VariableNames;
    ds.SelectedVariableNames = allVars([1 3 8 9]); % MNX_ID, reference, InChIKey, SMILES
    sel = ds.SelectedVariableNames;

    mnxMap   = containers.Map('KeyType','char', 'ValueType','any');
    chebiMap = containers.Map('KeyType','char', 'ValueType','any');

    while hasdata(ds)
        T = read(ds); if isempty(T), continue; end
        idCol = T.(sel{1}); refCol= T.(sel{2}); ikCol = T.(sel{3}); smCol = T.(sel{4});
        ok = (strlength(idCol) > 0) & ( (strlength(smCol) > 0) | (strlength(ikCol) > 0) );
        if ~any(ok), continue; end
        idCol = idCol(ok); refCol = refCol(ok); ikCol = ikCol(ok); smCol = smCol(ok);

        for k = 1:numel(idCol)
            key = char(idCol(k));
            if ~isKey(mnxMap, key)
                rec = struct('smiles', char(orEmpty(smCol(k))), ...
                             'inchikey', clean_inchikey(char(orEmpty(ikCol(k)))));
                mnxMap(key) = rec;
            end
        end
        for k = 1:numel(refCol)
            r = lower(char(refCol(k))); if isempty(r), continue; end
            toks = regexp(r, 'chebi:\d+', 'match'); if isempty(toks), continue; end
            rec = struct('smiles', char(orEmpty(smCol(k))), ...
                         'inchikey', clean_inchikey(char(orEmpty(ikCol(k)))));
            for t = 1:numel(toks)
                ck = toks{t}; if ~isKey(chebiMap, ck), chebiMap(ck) = rec; end
            end
        end
    end
end

function ids = extract_mnx_ids_from_miriam(miriam)
% Extract 'MNXM\d+' from miriam (prioritizing 'metanetx.chemical' entries)
    ids = {};
    mNames = field_to_cellstr(miriam, 'name');
    mVals  = field_to_cellstr(miriam, 'value');
    pri = find(contains(lower(mNames), 'metanetx.chemical'));
    for p = pri(:).'
        tok = regexp(mVals{p}, 'MNXM\d+', 'match');
        if ~isempty(tok), ids = [ids, tok]; end %#ok<AGROW>
    end
    if isempty(ids)
        for j = 1:numel(mVals)
            tok = regexp(mVals{j}, 'MNXM\d+', 'match');
            if ~isempty(tok), ids = [ids, tok]; end %#ok<AGROW>
        end
    end
    if ~isempty(ids), ids = unique(ids, 'stable'); end
end

function ids = extract_chebi_ids_from_miriam(miriam, toLower)
% Return as 'chebi:12345' if toLower=true; otherwise 'CHEBI:12345'
    ids = {};
    mNames = field_to_cellstr(miriam, 'name');
    mVals  = field_to_cellstr(miriam, 'value');
    emit = @(num) iff(toLower, ['chebi:' num], ['CHEBI:' num]);

    pri = find(contains(lower(mNames), 'chebi'));
    for p = pri(:).'
        s = mVals{p};
        tok = regexp(s, '(?i)CHEBI:\s*(\d+)', 'tokens');
        for t = 1:numel(tok), ids{end+1} = emit(tok{t}{1}); end %#ok<AGROW>
        tok2 = regexp(s, '^\s*(\d+)\s*$', 'tokens', 'once');
        if ~isempty(tok2), ids{end+1} = emit(tok2{1}); end %#ok<AGROW>
    end
    for j = 1:numel(mVals)
        tok = regexp(mVals{j}, '(?i)CHEBI:\s*(\d+)', 'tokens');
        for t = 1:numel(tok), ids{end+1} = emit(tok{t}{1}); end %#ok<AGROW>
    end
    if ~isempty(ids), ids = unique(ids, 'stable'); end
end

function changed = upsert_metinfo(dbMap, name, smiles, inchikey, verbosity)
% Upsert (name -> {smiles, inchikey}) into dbMap; fill-missing-only; normalize InChIKey.
    changed  = false;
    smiles   = orEmpty(smiles);
    inchikey = clean_inchikey(orEmpty(inchikey));
    if isempty(name), return; end

    if isKey(dbMap, name)
        rec = dbMap(name);
        oldS = orEmpty(rec.smiles);
        oldI = clean_inchikey(orEmpty(rec.inchikey));

        newS = oldS; newI = oldI;
        if isempty(oldS) && ~isempty(smiles),   newS = smiles;   end
        if isempty(oldI) && ~isempty(inchikey), newI = inchikey; end

        if ~strcmp(orEmpty(rec.smiles), newS) || ~strcmp(clean_inchikey(orEmpty(rec.inchikey)), newI)
            dbMap(name) = struct('smiles', newS, 'inchikey', newI);
            changed = true;
        end
    else
        if ~isempty(smiles) || ~isempty(inchikey)
            dbMap(name) = struct('smiles', smiles, 'inchikey', inchikey);
            changed = true;
        end
    end
end

function save_metinfo(dbMap, metInfoFile)
% Atomically rewrite metInfo.tsv from dbMap (unique names, no duplicates).
    [parentDir,~,~] = fileparts(metInfoFile);
    if ~isempty(parentDir) && ~exist(parentDir, 'dir'), mkdir(parentDir); end
    K = dbMap.keys; K = sort(K);

    tmpFile = [metInfoFile '.tmp'];
    fid = fopen(tmpFile, 'w');
    if fid < 0, error('Failed to open temp file for writing: %s', tmpFile); end

    for i = 1:numel(K)
        nm = K{i}; rec = dbMap(nm);
        sm  = orEmpty(rec.smiles);
        ik  = clean_inchikey(orEmpty(rec.inchikey));
        fprintf(fid, '%s\t%s\t%s\n', nm, sm, ik);
    end

    st = fclose(fid);
    if st ~= 0, error('Failed to close temp file: %s', tmpFile); end
    if exist(tmpFile, 'file') ~= 2, error('Temp file missing after write: %s', tmpFile); end

    if exist(metInfoFile, 'file') == 2
        try delete(metInfoFile); catch, end
    end
    [ok, msg] = movefile(tmpFile, metInfoFile, 'f');
    if ~ok
        [ok2, msg2] = copyfile(tmpFile, metInfoFile, 'f');
        if ok2, try delete(tmpFile); end %#ok<TRYNC>
        else, error('Failed to replace metInfo.tsv. movefile: %s | copyfile: %s', msg, msg2);
        end
    end
end

% -------- tiny utilities --------
function v = iff(cond, a, b), if cond, v = a; else, v = b; end, end
function v = ifthen(cond, a, b), if cond, v = a; else, v = b; end, end

function s = orEmpty(x)
    try
        if isempty(x), s = ''; return; end
        if iscell(x)
            if isempty(x{1}), s = ''; else, s = orEmpty(x{1}); end
            return;
        end
        if isstring(x)
            x = x(:);
            mask = ~ismissing(x) & (strlength(x) > 0);
            if any(mask), s = char(x(find(mask,1,'first'))); else, s = ''; end
            return;
        end
        if ischar(x), s = x; return; end
        s = char(string(x));
        if all(strlength(string(s)) == 0), s = ''; end
    catch
        s = '';
    end
end

function arr = field_to_cellstr(s, fname)
    arr = {};
    if ~isfield(s, fname), return; end
    v = s.(fname);
    if iscell(v)
        arr = cellfun(@char, v, 'UniformOutput', false);
    else
        arr = {char(v)};
    end
end

function t0 = rate_limit(t0, minInterval)
    if isempty(t0), t0 = tic; return; end
    dt = toc(t0);
    if dt < minInterval, pause(minInterval - dt); end
    t0 = tic;
end

function logv(verbosity, lvl, fmt, varargin)
    if verbosity >= lvl, fprintf(fmt, varargin{:}); end
end

function s = safeShort(s)
    MAXL = 60;
    if isempty(s), return; end
    if numel(s) > MAXL, s = [s(1:MAXL) '...']; end
end

% -------- ChEBI helpers (STEP B) --------
function enc = encode_chebi_for_path(idstr)
    enc = strrep(string(idstr), ":", "%3A");
    enc = regexprep(enc, ' ', '%20');
end

function J = safe_jsondecode(txt)
    J = [];
    try J = jsondecode(char(txt)); catch, end
end

function [smiles, inchikey] = extract_smiles_inchikey(J, rawTxt)
% Walk struct; fallback to regex. InChIKey normalization is done by caller.
    smiles = ""; inchikey = "";
    if ~isempty(J) && isstruct(J)
        [smiles, inchikey] = walk_for_keys(J);
        if smiles~="" || inchikey~="", return; end
    end
    s = char(string(rawTxt));
    m1 = regexp(s,'"(?:canonical_?smiles|smiles)"\s*:\s*"([^"]+)"','tokens','once');
    if ~isempty(m1), smiles = string(m1{1}); end
    m2 = regexp(s,'"(?:standard_?inchi_?key|inchi_?key)"\s*:\s*"([^"]+)"','tokens','once');
    if ~isempty(m2), inchikey = string(m2{1}); end
end

function [sm, ik] = walk_for_keys(S)
    sm = ""; ik = "";
    if ~isstruct(S), return; end
    fns = fieldnames(S);
    for i=1:numel(fns)
        fn = fns{i}; val = S.(fn); lfn = lower(fn);
        if ischar(val) || isstring(val)
            if contains(lfn,'smiles') && sm=="", sm = string(val); end
            if contains(lfn,'inchi') && contains(lfn,'key') && ik=="", ik = string(val); end
        elseif isstruct(val)
            [sm2, ik2] = walk_for_keys(val);
            if sm=="", sm = sm2; end
            if ik=="", ik = ik2; end
        elseif iscell(val)
            for k=1:numel(val)
                [sm2, ik2] = walk_for_keys(val{k});
                if sm=="", sm = sm2; end
                if ik=="", ik = ik2; end
            end
        end
        if sm~="" && ik~="", return; end
    end
end

function y = clean_inchikey(x)
% Strip optional 'InChIKey=' prefix (case-insensitive), trim spaces.
% Example input: 'InChIKey=BSYNRYMUTXBXSQ-UHFFFAOYSA-N' -> 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
    y = orEmpty(x);
    if isempty(y), return; end
    y = regexprep(y, '^\s*(?i)InChIKey\s*=\s*', '');
    y = strtrim(y);
end
