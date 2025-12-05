function mapTbl = map_basic_one_enzyme_per_row(model_basic, model_integrated)
% MAP_BASIC_ONE_ENZYME_PER_ROW
% Map BASIC vs INTEGRATED enzyme-constraint rows by:
%   key = <baseRxnID> + "|" + <enzyme multiset (order-insensitive, multiplicity-sensitive)>
%
% Multiset construction:
%   - For each nonzero enzyme column in rxnEnzMat(r,:), collect its UniProt token(s)
%   - Each token is uppercased and trimmed
%   - If the coefficient is an integer > 1 (within tolerance), replicate the token that many times
%   - The multiset is sorted to make a canonical signature (order-insensitive, multiplicity kept)
%
% BASIC normalization  : remove exactly a leading '^\d{3}_' (001_/002_/003_/...)
% INTEGRATED normalization: remove exactly a trailing '_EXP_\d+$'
%
% Output columns:
%   baseID                 normalized reaction ID used for matching
%   enzymeBag_basic        1xK cellstr (tokens, sorted; multiplicity preserved)
%   enzymeSig_basic        canonical signature string of the bag (e.g. 'P00561|P00561')
%   basic_rxnID            ec rxnID in BASIC enzymeConstraints
%   basic_rxnName          resolved name from model_basic.rxns (if available)
%   basic_row              row index in basic enzymeConstraints.rxnEnzMat
%
%   enzymeBag_integrated   1xK cellstr for the matched INTEGRATED row (first match)
%   enzymeSig_integrated   canonical signature string for the matched INTEGRATED row
%   integrated_rxnID       matched ec rxnID in INTEGRATED enzymeConstraints (first match)
%   integrated_rxnName     resolved name (if available)
%   integrated_row         row index in integrated enzymeConstraints.rxnEnzMat (first match)
%
%   match_count            number of INTEGRATED rows sharing the same (baseID, bag) key
%   integrated_candidates  cell array of all candidate INTEGRATED rxnIDs (if >1)
%   status                 'matched' | 'no_match'
%   note                   diagnostic text
%
% Assumptions:
%   - model_*.enzymeConstraints.{rxns,enzymes,rxnEnzMat} exist
%   - enzymeConstraints.enzymes entries may be char/string or cells of tokens
%   - rxnEnzMat can be sparse; coefficients may be non-integer
%
% Author: ChatGPT (per your exact matching spec)

    % -------- sanity checks --------
    mustHave = {'enzymeConstraints'};
    for i = 1:numel(mustHave)
        if ~isfield(model_basic, mustHave{i}) || ~isfield(model_integrated, mustHave{i})
            error('Both inputs must contain .enzymeConstraints.');
        end
    end
    reqEC = {'rxns','enzymes','rxnEnzMat'};
    for i = 1:numel(reqEC)
        if ~isfield(model_basic.enzymeConstraints, reqEC{i})
            error('model_basic.enzymeConstraints.%s is required.', reqEC{i});
        end
        if ~isfield(model_integrated.enzymeConstraints, reqEC{i})
            error('model_integrated.enzymeConstraints.%s is required.', reqEC{i});
        end
    end

    bEC = model_basic.enzymeConstraints;
    iEC = model_integrated.enzymeConstraints;

    b_rxn = string(bEC.rxns(:));
    i_rxn = string(iEC.rxns(:));

    % -------- normalize reaction IDs for matching --------
    baseB = stripBasicPrefix(b_rxn);        % remove '001_'
    baseI = stripIntegratedSuffix(i_rxn);   % remove '_EXP_1'

    % -------- build enzyme multisets & canonical signatures --------
    [bagB, sigB] = enzyme_bag_and_signature(bEC.rxnEnzMat, bEC.enzymes);
    [bagI, sigI] = enzyme_bag_and_signature(iEC.rxnEnzMat, iEC.enzymes);

    % -------- build index for INTEGRATED: (baseID, signature) → row(s) --------
    keyI_map = containers.Map('KeyType','char','ValueType','any');
    for r = 1:numel(i_rxn)
        if strlength(sigI{r}) == 0 || strlength(baseI(r)) == 0
            % skip rows without any enzyme tokens or empty normalized base
            continue;
        end
        k = char(baseI(r) + "|" + string(sigI{r}));
        if isKey(keyI_map, k)
            keyI_map(k) = [keyI_map(k), r]; %#ok<AGROW>
        else
            keyI_map(k) = r;
        end
    end

    % -------- map BASIC → INTEGRATED --------
    nB = numel(b_rxn);

    baseID               = baseB;
    enzymeBag_basic      = bagB;      % 1xK cell arrays (per row)
    enzymeSig_basic      = strings(nB,1);
    basic_rxnID          = b_rxn;
    basic_rxnName        = resolve_names_from_model(model_basic, b_rxn);
    basic_row            = (1:nB).';

    enzymeBag_integrated = repmat({cell(1,0)}, nB, 1);
    enzymeSig_integrated = strings(nB,1);
    integrated_rxnID     = strings(nB,1);
    integrated_rxnName   = strings(nB,1);
    integrated_row       = nan(nB,1);

    match_count          = zeros(nB,1);
    integrated_candidates= cell(nB,1);
    status               = strings(nB,1);
    note                 = strings(nB,1);

    for r = 1:nB
        enzymeSig_basic(r) = string(sigB{r});
        if strlength(enzymeSig_basic(r)) == 0
            status(r) = "no_match";
            note(r)   = "basic row has empty enzyme multiset";
            continue;
        end
        if strlength(baseID(r)) == 0
            status(r) = "no_match";
            note(r)   = "basic row has empty normalized base reaction ID";
            continue;
        end

        key = char(baseID(r) + "|" + enzymeSig_basic(r));
        if ~isKey(keyI_map, key)
            status(r) = "no_match";
            note(r)   = "no integrated row with same <baseID|enzyme multiset>";
            continue;
        end

        cand = keyI_map(key);
        match_count(r) = numel(cand);
        integrated_candidates{r} = cellstr(i_rxn(cand));

        % Choose the first candidate as canonical; keep others in candidates
        ii = cand(1);
        integrated_row(r)     = ii;
        integrated_rxnID(r)   = i_rxn(ii);
        integrated_rxnName(r) = resolve_name_from_model(model_integrated, i_rxn(ii));
        enzymeBag_integrated{r} = bagI{ii};
        enzymeSig_integrated(r) = string(sigI{ii});

        status(r) = "matched";
        if numel(cand) > 1
            note(r) = sprintf('multiple integrated rows share key (%d candidates)', numel(cand));
        else
            note(r) = "";
        end
    end

    % -------- assemble output table --------
    mapTbl = table( ...
        baseID, ...
        enzymeBag_basic, enzymeSig_basic, ...
        basic_rxnID, basic_rxnName, basic_row, ...
        enzymeBag_integrated, enzymeSig_integrated, ...
        integrated_rxnID, integrated_rxnName, integrated_row, ...
        match_count, integrated_candidates, ...
        status, note ...
    );
end

% ====================== helpers ======================

function s = stripBasicPrefix(ids)
% BASIC: remove exactly a leading '^\d{3}_' (001_/002_/003_/...)
    s = string(ids);
    s = regexprep(s, '^[0-9]{3}_', '');
end

function s = stripIntegratedSuffix(ids)
% INTEGRATED: remove exactly a trailing '_EXP_<n>'
    s = string(ids);
    s = regexprep(s, '_EXP_[0-9]+$', '');
end

function names = resolve_names_from_model(model, rxnIDs)
% Vectorized resolution of rxn names from model.rxns → model.rxnNames
    names = strings(numel(rxnIDs),1);
    if ~isfield(model,'rxns') || ~isfield(model,'rxnNames'), return; end
    mIDs = string(model.rxns(:));
    mNm  = string(model.rxnNames(:));
    [ok, pos] = ismember(string(rxnIDs(:)), mIDs);
    names(ok) = mNm(pos(ok));
end

function name = resolve_name_from_model(model, rxnID)
% Scalar convenience version
    name = "";
    if ~isfield(model,'rxns') || ~isfield(model,'rxnNames'), return; end
    idx = find(strcmp(model.rxns, rxnID), 1, 'first');
    if ~isempty(idx), name = string(model.rxnNames{idx}); end
end

function [bags, sigs] = enzyme_bag_and_signature(rxnEnzMat, enzymes)
% Build per-row enzyme bags (1xK cellstr) and canonical signatures (char)
% - tokens are uppercase UniProt IDs
% - multiplicity preserved; integer stoichiometry >1 → replicate that many times
% - signature is a '|'-joined string of sorted tokens (e.g., 'P00561|Q9XYZ1|Q9XYZ1')

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

                if abs(coeff) <= tol
                    continue; % zero within tolerance
                end

                % Pull raw enzyme entry; may be char/string or a cell of tokens
                raw = enzymes{c};
                tokens = normalize_enzyme_entry(raw);

                if isempty(tokens)
                    % skip silently (no valid token)
                    continue;
                end

                % multiplicity: replicate if coefficient is a positive integer
                rep = 1;
                if coeff > 0
                    zr = abs(coeff - round(coeff));
                    if zr <= 1e-6 && round(coeff) >= 1
                        rep = round(coeff);
                    end
                end

                % append tokens with multiplicity
                for t = 1:numel(tokens)
                    toks = [toks, repmat({tokens{t}}, 1, rep)]; %#ok<AGROW>
                end
            end
        end

        % Sort to build canonical, order-insensitive signature
        if isempty(toks)
            bags{r} = cell(1,0);
            sigs{r} = '';
        else
            toks = sort(toks);
            bags{r} = toks;
            sigs{r} = strjoin(toks, '|'); % e.g., 'P00561|P00561'
        end
    end
end

function tokens = normalize_enzyme_entry(raw)
% Convert enzyme entry to a flat list of uppercase UniProt tokens.
% Accepts:
%   - char/string: one token
%   - cell: list of tokens (each char/string)
% Drops empties. Does NOT split on whitespace inside a token.

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
