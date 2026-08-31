function model = getECnumber(model, parameters)
% getECnumber  Populate/validate EC numbers for enzyme-constrained model
%   1) Pull ECs from model (map from model.eccodes)
%   2) Backfill missing/invalid ECs from local UniProt (uniprot.tsv parsed by ParseUniProtData)
%
% Final ECs are written to model.enzymeConstraints.eccodes (aligned to enzymeConstraints.rxns).

    % -------- Args / guards --------
    if nargin < 2 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end
    if ~isfield(model, 'enzymeConstraints') || isempty(model.enzymeConstraints)
        error('The model does not contain the enzymeConstraints structure.');
    end

    EC = model.enzymeConstraints;
    needFields = {'rxns','enzymes','rxnEnzMat','ecModeltype'};
    for f = needFields
        if ~isfield(EC, f{1})
            error('model.enzymeConstraints.%s is missing.', f{1});
        end
    end

    % -------- Load UniProt (local cache or download) --------
    uniprot_Path = fullfile(parameters.reconstructionDir, 'uniprot.tsv');
    if ~isfile(uniprot_Path)
        DownloadUniProtData(parameters.uniprot, parameters.reconstructionDir);
    end
    dbStruct = ParseUniProtData(uniprot_Path);
    if isempty(dbStruct) || ~isfield(dbStruct,'ID') || ~isfield(dbStruct,'eccodes')
        error('Parsed UniProt TSV missing required fields {ID, eccodes}.');
    end

    % -------- Reaction name normalization --------
    switch lower(EC.ecModeltype)
        case {'integrated','isozyme'}
            rxnNames = EC.rxns;
        otherwise
            rxnNames = extractAfter(EC.rxns, 4); % strip 4-char prefix like 'prot'
    end

    % Map enzymeConstraints.rxns to model.rxns only when needed (never index tf by idx!)
    [rxnFoundInModel, rxnPosInModel] = ismember(rxnNames, model.rxns);

    % -------- Prepare "source EC" from model --------
    nRxn = numel(EC.rxns);
    outEC = repmat({''}, nRxn, 1);     % final ECs aligned to enzymeConstraints.rxns

    % Prefer enzymeConstraints.eccodes if shape matches
    if isfield(EC, 'eccodes') && numel(EC.eccodes) == nRxn
        srcEC = EC.eccodes;
    else
        % Otherwise, try map from model.eccodes (length must equal model.rxns)
        if isfield(model, 'eccodes') && ~isempty(model.eccodes)
            srcEC_all = asCellstr(model.eccodes);
            if numel(srcEC_all) ~= numel(model.rxns)
                warning('model.eccodes length (%d) != model.rxns length (%d); ignoring model.eccodes.', ...
                        numel(srcEC_all), numel(model.rxns));
                srcEC = repmat({''}, nRxn, 1);
            else
                srcEC = repmat({''}, nRxn, 1);
                hit = rxnFoundInModel;                 % only where the rxn maps
                pos = rxnPosInModel(hit);              % positions in model.rxns
                srcEC(hit) = srcEC_all(pos);           % map onto enzymeConstraints order
            end
        else
            srcEC = repmat({''}, nRxn, 1);
        end
    end
%%
    % Normalize source EC strings
    srcEC = normalizeECStrings(srcEC);

    % -------- Validate EC format (full-string match; allow multiple ';') --------
    validPattern = '^(?:\d+(?:\.(?:\d+|-)){3})(?:;(?:\d+(?:\.(?:\d+|-)){3}))*$';
    isValidSrc = ~cellfun(@isempty, srcEC) & ~cellfun(@isempty, regexp(srcEC, validPattern, 'once'));
    outEC(isValidSrc) = srcEC(isValidSrc);

    % -------- Backfill from UniProt for missing/invalid --------
    % Prepare UniProt ID -> EC mapping
    dbIDs   = dbStruct.ID(:);
    dbECraw = normalizeECStrings(dbStruct.eccodes(:));

    enzList = EC.enzymes(:);  % UniProt accessions used in the model
    [idFound, idLoc] = ismember(enzList, dbIDs);
    enz2ec = repmat({''}, numel(enzList), 1);
    enz2ec(idFound) = dbECraw(idLoc(idFound));

    % For reactions still lacking EC, collect from enzymes (columns with nonzero rxnEnzMat)
    needFill = find(cellfun(@isempty, outEC));
    for i = needFill(:)'
        cols = find(EC.rxnEnzMat(i, :) ~= 0);
        if isempty(cols), continue; end
        cands = enz2ec(cols);
        cands = cands(~cellfun(@isempty, cands));
        if isempty(cands), continue; end

        % Split multi-EC strings, deduplicate, validate each token
        toks = cellfun(@(s) strsplit(s, ';'), cands, 'UniformOutput', false);
        toks = [toks{:}];
        toks = unique(strtrim(toks));
        toks = toks(~cellfun(@isempty, toks));
        if isempty(toks), continue; end

        good = ~cellfun(@isempty, regexp(toks, '^\d+(?:\.(?:\d+|-)){3}$', 'once'));
        toks = toks(good);
        if isempty(toks), continue; end

        outEC{i} = strjoin(toks, ';');
    end

    % -------- Write back --------
    model.enzymeConstraints.eccodes = outEC;

    % -------- Report --------
    nAssigned = sum(~cellfun(@isempty, outEC));
    nFromModel = sum(isValidSrc);
    fprintf('[getECnumber] ECs assigned for %d/%d reactions (from model: %d, from UniProt: %d).\n', ...
        nAssigned, nRxn, nFromModel, nAssigned - nFromModel);
end

% ==== helpers ====

function ec = asCellstr(x)
% convert char/string/cellstr to column cellstr
    if isstring(x), x = cellstr(x); end
    if ischar(x),   x = cellstr(x); end
    if ~iscell(x),  x = {x}; end
    x = x(:);
    ec = cell(size(x));
    for k = 1:numel(x)
        if isstring(x{k}), ec{k} = char(x{k});
        elseif ischar(x{k}), ec{k} = x{k};
        else, ec{k} = char(string(x{k}));
        end
    end
end

function ec = normalizeECStrings(ec)
% Normalize EC strings:
%   - convert to cellstr (column)
%   - trim spaces
%   - remove optional leading 'EC'
%   - normalize spaces around '.' and ';'
    ec = asCellstr(ec);
    ec = cellfun(@strtrim, ec, 'UniformOutput', false);
    ec = regexprep(ec, '^\s*EC\s*', '', 'ignorecase');
    ec = regexprep(ec, '\s*\.\s*', '.');
    ec = regexprep(ec, '\s*;\s*', ';');
end
