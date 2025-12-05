function [nameMap, matchTbl] = map_ec_rxnname_intersect(models)
% MAP_EC_RXNNAME_INTERSECT
% Intersect ecModel enzyme-constraint rows from 2 or 3 models using:
%   key = <normalized base rxnID> + "|" + <enzyme multiset signature>
%
% Normalization by ecModeltype (auto-detected from model.enzymeConstraints.ecModeltype):
%   - 'basic'     : strip exactly a leading 001_/002_/...  (regex ^\d{3}_)
%   - 'integrated': strip exactly a trailing _EXP_<n>      (regex _EXP_\d+$)
%   - 'isozyme'   : SAME as 'integrated' (strip _EXP_<n>)
%
% Enzyme multiset signature (order-insensitive, multiplicity-sensitive):
%   - Built from rxnEnzMat row using EC.enzymes (UniProt IDs)
%   - Uppercased tokens, integer stoich >1 is replicated
%   - Signature = sorted tokens joined by '|'
%
% Inputs
%   models : cell array of 2 or 3 ecModels (structs with enzymeConstraints)
%
% Outputs
%   nameMap : table with ONLY matched EC-layer rxnIDs (2 or 3 columns),
%             column names are derived from ecModeltype: 'basic','integrated','isozyme'
%   matchTbl: table summarizing matches per (baseID, enzymeSig) key.
%
% Notes
% - Fixes your error by storing labels via cell braces {} and sanitizing names.
% - Dynamic field and table var names use char for broader MATLAB compatibility.

    % -------- arguments & checks --------
    if ~iscell(models) || numel(models) < 2 || numel(models) > 3
        error('models must be a cell array with 2 or 3 ecModels.');
    end

    nM    = numel(models);
    M     = cell(1,nM);          % per-model cache
    keysAll = cell(1,nM);        % all keys per model (char cell)
    labels  = cell(1,nM);        % column labels (sanitized, char)

    % -------- preprocess each model --------
    for m = 1:nM
        mdl = models{m};
        ensure_ec_fields(mdl, m);

        EC  = mdl.enzymeConstraints;
        rxn = string(EC.rxns(:));

        % detect ecModeltype and sanitize label
        if isfield(EC,'ecModeltype') && ~isempty(EC.ecModeltype)
            kind = lower(char(string(EC.ecModeltype)));
        end
        % store label into CELL using braces; make it a valid MATLAB name (char)
        labels{m} = char(matlab.lang.makeValidName(kind));

        % normalize rxnIDs by kind
        switch kind
            case 'basic'
                baseID = strip_basic_prefix(rxn);
            case {'integrated','isozyme'}
                baseID = strip_integrated_suffix(rxn);
            otherwise
                % unknown tag -> safest is integrated-style suffix strip
                baseID = strip_integrated_suffix(rxn);
        end

        % enzyme multiset & signatures
        [bag, sig] = enzyme_bag_and_signature(EC.rxnEnzMat, EC.enzymes);

        % build key->rows map (skip empties)
        keyMap = containers.Map('KeyType','char','ValueType','any');
        modelKeys = cell(0,1);
        for r = 1:numel(rxn)
            if strlength(baseID(r))==0 || isempty(sig{r}) || strlength(sig{r})==0
                continue;
            end
            k = char(baseID(r) + "|" + string(sig{r}));
            modelKeys{end+1,1} = k; %#ok<AGROW>
            if isKey(keyMap,k)
                keyMap(k) = [keyMap(k), r]; %#ok<AGROW>
            else
                keyMap(k) = r;
            end
        end

        % cache
        M{m}.rxn    = rxn;
        M{m}.baseID = baseID;
        M{m}.bag    = bag;
        M{m}.sig    = sig;
        M{m}.keyMap = keyMap;
        keysAll{m}  = unique(modelKeys);
    end

    % -------- compute intersection of keys across provided models --------
    commonKeys = keysAll{1};
    for m = 2:nM
        commonKeys = intersect(commonKeys, keysAll{m});
    end
    commonKeys = sort(commonKeys);

    % -------- assemble nameMap (only matched names) --------
    % Choose the FIRST candidate row per model for each key (stable).
    nmCols = cell(1,nM);
    nmCols(:) = {strings(numel(commonKeys),1)};
    for i = 1:numel(commonKeys)
        k = commonKeys{i};
        for m = 1:nM
            rows = M{m}.keyMap(k);
            nmCols{m}(i) = M{m}.rxn(rows(1));
        end
    end

    % Build the table with exactly 2 or 3 name columns
    Tnm = table();
    for m = 1:nM
        % labels{m} is char; safe as a table variable name
        Tnm.(labels{m}) = nmCols{m};
    end
    nameMap = Tnm;

    % -------- assemble matchTbl (summary per match key) --------
    % Split key -> [baseID, enzymeSig]
    baseID_vec  = strings(numel(commonKeys),1);
    enzSig_vec  = strings(numel(commonKeys),1);
    for i = 1:numel(commonKeys)
        [b, s] = split_key(commonKeys{i});
        baseID_vec(i) = b;
        enzSig_vec(i) = s;
    end

    % Per-model chosen row, chosen rxnID, candidate counts & lists
    chosen_row  = cell(1,nM);
    chosen_rxn  = cell(1,nM);
    cand_count  = cell(1,nM);
    cand_list   = cell(1,nM);
    for m = 1:nM
        chosen_row{m} = nan(numel(commonKeys),1);
        chosen_rxn{m} = strings(numel(commonKeys),1);
        cand_count{m} = zeros(numel(commonKeys),1);
        cand_list{m}  = cell(numel(commonKeys),1);
    end

    for i = 1:numel(commonKeys)
        k = commonKeys{i};
        for m = 1:nM
            rows = M{m}.keyMap(k);
            cand_count{m}(i) = numel(rows);
            cand_list{m}{i}  = cellstr(M{m}.rxn(rows));
            chosen_row{m}(i) = rows(1);
            chosen_rxn{m}(i) = M{m}.rxn(rows(1));
        end
    end

    % Build matchTbl with dynamic per-model columns
    matchTbl = table(baseID_vec, enzSig_vec, 'VariableNames', {'baseID','enzymeSig'});
    for m = 1:nM
        % use labels{m} (char) to build dynamic field names
        matchTbl.([labels{m} '_rxnID'])       = chosen_rxn{m};
        matchTbl.([labels{m} '_row'])         = chosen_row{m};
        matchTbl.([labels{m} '_candCount'])   = cand_count{m};
        matchTbl.([labels{m} '_candList'])    = cand_list{m};
    end
end

% ====================== local helpers ======================

function ensure_ec_fields(model, mpos)
    if ~isfield(model,'enzymeConstraints')
        error('models{%d} is missing .enzymeConstraints.', mpos);
    end
    EC = model.enzymeConstraints;
    need = {'rxns','enzymes','rxnEnzMat'};
    for i = 1:numel(need)
        if ~isfield(EC, need{i})
            error('models{%d}.enzymeConstraints.%s is required.', mpos, need{i});
        end
    end
end

function s = strip_basic_prefix(ids)
% Remove exactly a leading '001_/002_/...' (^\d{3}_)
    s = string(ids);
    s = regexprep(s, '^[0-9]{3}_', '');
end

function s = strip_integrated_suffix(ids)
% Remove exactly a trailing '_EXP_<n>'
    s = string(ids);
    s = regexprep(s, '_EXP_[0-9]+$', '');
end

function [bags, sigs] = enzyme_bag_and_signature(rxnEnzMat, enzymes)
% Build per-row enzyme bags (1xK cellstr) and canonical signatures (char).
% - tokens are uppercase UniProt IDs
% - multiplicity preserved; integer stoichiometry >1 → replicate
% - signature = sorted tokens joined with '|'
    tol = 1e-9;
    nR  = size(rxnEnzMat,1);
    [rIdx, cIdx, v] = find(rxnEnzMat);

    rows2cols = accumarray(rIdx, cIdx, [nR,1], @(x){x}, {[]});
    rows2vals = accumarray(rIdx, v,    [nR,1], @(x){x}, {[]});

    bags = cell(nR,1);
    sigs = cell(nR,1);

    for r = 1:nR
        cols = rows2cols{r};
        vals = rows2vals{r};
        toks = {};
        if ~isempty(cols)
            for k = 1:numel(cols)
                c = cols(k);
                coeff = vals(k);
                if abs(coeff) <= tol, continue; end

                raw = enzymes{c};                     % may be char/string or cell
                tokens = normalize_enzyme_entry(raw); % -> cellstr uppercase
                if isempty(tokens), continue; end

                rep = 1;
                if coeff > 0
                    zr = abs(coeff - round(coeff));
                    if zr <= 1e-6 && round(coeff) >= 1
                        rep = round(coeff);
                    end
                end
                for t = 1:numel(tokens)
                    toks = [toks, repmat({tokens{t}}, 1, rep)]; %#ok<AGROW>
                end
            end
        end

        if isempty(toks)
            bags{r} = cell(1,0);
            sigs{r} = '';
        else
            toks = sort(toks);          % order-insensitive
            bags{r} = toks;
            sigs{r} = strjoin(toks, '|');
        end
    end
end

function tokens = normalize_enzyme_entry(raw)
% Convert enzyme entry to a flat list of uppercase UniProt tokens.
% Accepts char/string or cell of those; drops empties.
    if iscell(raw)
        flat = raw(:);
    else
        flat = {raw};
    end
    tokens = {};
    for i = 1:numel(flat)
        try
            s = strtrim(char(string(flat{i})));
        catch
            s = '';
        end
        if isempty(s), continue; end
        tokens{end+1} = upper(s); %#ok<AGROW>
    end
end

function [b, s] = split_key(k)
% Split "baseID|signature" into its parts. Signature may contain more '|'s.
    parts = regexp(k, '\|', 'split');
    if numel(parts) >= 2
        b = string(parts{1});
        s = string(strjoin(parts(2:end), '|'));
    else
        b = string(k);
        s = "";
    end
end
