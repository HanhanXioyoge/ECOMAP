function writtenTable = writeInputFile(model, DeepLearningModel, parameters)
% WRITEINPUTFILE
% Build input CSV files for selected deep-learning kcat models (DLKcat, UniKP, CatPred)
% from an enzyme-constrained model (ecModel).
%
% The function enumerates (substrate, reaction, protein) tuples from the
% enzymeConstraints mapping, filters out ignored/currency metabolites, and writes
% per-model CSVs with the columns each predictor typically expects.
%
% Inputs
%   model              : COBRA/GECKO-like model with model.enzymeConstraints.* fields
%   DeepLearningModel  : char/string/array/cellstr, subset of {'DLKcat','UniKP','CatPred'}
%   parameters         : struct with fields
%       - reconstructionDir          : output folder for CSVs (required)
%       - onlyWithSmiles   : logical (default true) keep only rows with SMILES
%       - org_name         : organism name to stamp in the table (default 'NA')
%
% Output
%   writtenTable       : the master table prior to per-model column selection
%
% Notes
%   - If model.metSmiles/metInChIKey/metMetaNetXID are missing, the function
%     degrades gracefully and fills 'NA' while warning.
%   - enzymeConstraints.PDB is optional; if absent, CatPred's pdbpath column is empty.

    % ------------------------- Parameter defaults -------------------------
    if nargin < 3 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.');
        end
    end
    if ~isfield(parameters, 'reconstructionDir') || isempty(parameters.reconstructionDir)
        error('parameters.reconstructionDir is required.');
    end
    reconstructionDir = parameters.reconstructionDir;

    if ~isfield(parameters, 'onlyWithSmiles') || isempty(parameters.onlyWithSmiles)
        parameters.onlyWithSmiles = true;
    end

    if ~isfield(parameters, 'org_name') || isempty(parameters.org_name)
        parameters.org_name = 'NA';
    end

    % ------------------------- Model list handling ------------------------
    if nargin < 2 || isempty(DeepLearningModel)
        error('writeInputFile:MissingDLModelList', ...
              'DeepLearningModel must be a non-empty list. Allowed: {DLKcat, UniKP, CatPred}.');
    end

    if ischar(DeepLearningModel) || isstring(DeepLearningModel)
        DeepLearningModel = cellstr(DeepLearningModel);
    elseif ~iscell(DeepLearningModel)
        error('DeepLearningModel must be a char/string/string array/cellstr.');
    end

    % Normalize to upper-case unique tags
    dlModels = upper(strtrim(string(DeepLearningModel)));
    dlModels = unique(cellstr(dlModels));  % remove duplicates
    validModels = {'DLKCAT','UNIKP','CATPRED'};
    bad = setdiff(dlModels, validModels);
    if ~isempty(bad)
        error('Unrecognized model types: %s. Allowed: {DLKcat, UniKP, CatPred}.', strjoin(bad, ', '));
    end

    % ------------------------- Basic ecModel checks -----------------------
    if ~isfield(model, 'enzymeConstraints') || isempty(model.enzymeConstraints)
        error('The model does not contain the enzymeConstraints structure.');
    end
    ecz = model.enzymeConstraints;
    reqFields = {'ecModeltype','rxns','genes','eccodes','sequence','rxnEnzMat'};
    tf = isfield(ecz, reqFields); 
    if ~all(tf)
        missing = reqFields(~tf);
        error('Missing model.enzymeConstraints field(s): %s.', strjoin(missing, ', '));
    end

    % Optional commonly used fields
    hasSmiles= isfield(model, 'metSmiles');
    commonly_fields = {'metSmiles', 'metInChIKey', 'metMetaNetXID'};
    tf = isfield(model, commonly_fields); 
    if ~all(tf)
        missing = commonly_fields(~tf);
        error('Missing model field(s): %s. Check if the functions addMetMetaNetXID and getMetinfo are executed', strjoin(missing, ', '));
    end
    model.metSmiles = string(model.metSmiles(:));

    % --------------------- Normalize reaction name set --------------------
    % Some ecModel variants store ec rxn IDs with a prefix (e.g. 'prot_');
    % Here we map those back to model.rxns for indexing S matrix.
    switch lower(ecz.ecModeltype)
        case {'integrated','isozyme'}
            ecRxnNames = ecz.rxns;                 % identical to model.rxns
        otherwise
            % strip a 4-char prefix such as 'prot' from ec rxn names
            ecRxnNames = extractAfter(ecz.rxns, 4);
    end


    [okMap, RxnIdxs] = ismember(ecRxnNames, model.rxns);
    if ~all(okMap)
        miss = ecRxnNames(~okMap);
        error('Not all enzymeConstraints.rxns could be mapped to model.rxns. Example missing: %s', miss{1});
    end

    % ------------------ Build a substrate-filtered S matrix ----------------
    % 1) Remove "ignored metabolites" (names and optional SMILES)
    metsNoSpecialChars = lower(regexprep(model.metNames, '[^0-9a-zA-Z]+', ''));
    fIgnore = fullfile(findECOMAProot, 'scripts', 'database', 'IgnoreMets.tsv');
    fID = fopen(fIgnore, 'r');
    if fID == -1
        error('Cannot open %s', fIgnore);
    end
    fileData = textscan(fID, '%s %s', 'delimiter', '\t');
    fclose(fID);
    ignoreMets   = lower(regexprep(fileData{1}, '[^0-9a-zA-Z]+', ''));
    ignoreSmiles = fileData{2};
    ignoreSmiles(cellfun(@isempty, ignoreSmiles)) = [];

    ignoreMetsIdx = ismember(metsNoSpecialChars, ignoreMets);
    if hasSmiles
        ignoreMetsIdx = ignoreMetsIdx | ismember(model.metSmiles, ignoreSmiles);
    end
    % 2) Remove protein-usage pseudo-mets
    if isfield(model, 'mets')
        ignoreMetsIdx = ignoreMetsIdx | startsWith(model.mets, 'prot_');
    end

    reducedS = model.S;
    reducedS(ignoreMetsIdx, :) = 0;

    % 3) Remove "currency pairs" (only when both appear as sub & prod)
    fCurrency = fullfile(findECOMAProot, 'scripts', 'database', 'CurrencyMets.tsv');
    fID = fopen(fCurrency, 'r');
    if fID == -1
        error('Cannot open %s', fCurrency);
    end
    curData = textscan(fID, '%s %s', 'delimiter', '\t');
    fclose(fID);
    curA = lower(regexprep(curData{1}, '[^0-9a-zA-Z]+', ''));
    curB = lower(regexprep(curData{2}, '[^0-9a-zA-Z]+', ''));
    for i = 1:numel(curA)
        subs = strcmp(curA{i}, metsNoSpecialChars);
        prod = strcmp(curB{i}, metsNoSpecialChars);
        if any(subs) && any(prod)
            [~, subsRxns] = find(reducedS(subs, :));
            [~, prodRxns] = find(reducedS(prod, :));
            pairRxns = intersect(subsRxns, prodRxns);
            if ~isempty(pairRxns)
                tmpS = reducedS;
                tmpS([find(subs); find(prod)], pairRxns) = 0;
                keep = any(tmpS(:, pairRxns) < 0, 1); % keep only if substrates remain
                reducedS([find(subs); find(prod)], pairRxns(keep)) = 0;
            end
        end
    end

    % 4) Keep only ec reactions of interest (columns in RxnIdxs)
    clearedS = reducedS(:, RxnIdxs);

    % ---------------- Enumerate (substrate, reaction, protein) -------------
    % Substrates are negative entries in S (:, rxn); collect all pairs
    [metIdx, rxnLocalIdx] = find(clearedS < 0);      % indices within clearedS

    % Pre-allocate dynamic cell arrays (growth-friendly in MATLAB R2021+)
    rows_rxn_name       = {};
    rows_organism       = {};
    rows_gene           = {};
    rows_protein        = {};
    rows_ecnum          = {};
    rows_metanetxid     = {};
    rows_sub_name       = {};
    rows_sub_smiles     = {};
    rows_inchi_key      = {};
    rows_seq            = {};
    rows_struct_id      = {};

    for k = 1:numel(rxnLocalIdx)
        rxnE = rxnLocalIdx(k);           % index in ecz.rxns / ecRxnNames
        pList = find(ecz.rxnEnzMat(rxnE, :));     % proteins linked to this reaction
        if isempty(pList)
            continue; % no enzymes assigned -> skip
        end

        subName     = safeGetCell(model.metNames, metIdx(k), '');
        metanetxid  = safeGetCell(model.metMetaNetXID, metIdx(k), '');
        inchi_key   = safeGetCell(model.metInChIKey, metIdx(k), '');
        subSmiles   = 'NA';
        if hasSmiles
            val = safeGetCell(model.metSmiles, metIdx(k), '');
            if ~isempty(val), subSmiles = val; end
        end

        for pp = 1:numel(pList)
            p = pList(pp);

            geneID    = safeGetCell(ecz.genes, p, 'NA');
            proteinID = safeGetCell(ecz.enzymes, p, 'NA');
            seq       = safeGetCell(ecz.sequence, p, '');
            ecnum     = safeGetCell(ecz.eccodes, rxnE, '');
            pdb       = safeGetCell(ecz.PDB, p, '');
            
            if ~isempty(seq) || ~strcmpi(subSmiles, 'NA')
                rows_rxn_name{end+1,1}       = safeGetCell(ecz.rxns, rxnE,  'NA');
                rows_organism{end+1,1}       = parameters.org_name;
                rows_gene{end+1,1}           = geneID;
                rows_protein{end+1,1}        = proteinID;
                rows_ecnum{end+1,1}          = ecnum;
                rows_metanetxid{end+1,1}     = metanetxid;
                rows_sub_name{end+1,1}       = subName;
                rows_sub_smiles{end+1,1}     = subSmiles;
                rows_inchi_key{end+1,1}      = inchi_key;
                rows_seq{end+1,1}            = seq;
                rows_struct_id{end+1,1}      = pdb;
            end
        end
    end

    % ------------------------- Build master table -------------------------
    T = table( ...
        rows_rxn_name, rows_organism, rows_gene,rows_protein,rows_ecnum, rows_metanetxid,...
        rows_sub_name,rows_sub_smiles, rows_inchi_key, rows_seq,rows_struct_id, ...
        'VariableNames', {'ReactionName', 'Organism', 'GeneID', 'ProteinID', ...
        'EC Number', 'MetaNetXID', 'Substrate', 'SMILES', 'InChIKey', 'sequence', 'pdbpath'} );

    % Drop rows without SMILES if requested
    if parameters.onlyWithSmiles
        keep = ~cellfun(@isempty, T.SMILES) & ~strcmpi(T.SMILES, 'NA');
        T = T(keep, :);
    end

    % De-duplicate identical rows (stable to preserve order)
    if ~isempty(T)
        T = unique(T, 'rows', 'stable');
    end

    % ---------------------- Write per-model CSV files ----------------------
    for i = 1:numel(dlModels)
        tag = dlModels{i};
        switch tag
            case 'DLKCAT'
                TL = T(:, {'ReactionName', 'Organism', 'GeneID', 'ProteinID', ...
                    'EC Number', 'MetaNetXID', 'Substrate', 'SMILES', 'InChIKey', 'sequence'});
                outFile = fullfile(reconstructionDir, 'kcatData', 'DLKcat_input.csv');
                writetable(TL, outFile);
                fprintf('DLKcat input written to: %s\n', outFile);

            case 'UNIKP'
                % If UniKP does not need organism, drop it; keep EC if available
                TU = T(:, {'ReactionName', 'Organism', 'GeneID', 'ProteinID', ...
                    'EC Number', 'MetaNetXID','Substrate', 'SMILES', 'InChIKey', 'sequence'});
                outFile = fullfile(reconstructionDir, 'kcatData', 'UniKP_input.csv');
                writetable(TU, outFile);
                fprintf('UniKP input written to: %s\n', outFile);

            case 'CATPRED'
                % CatPred requires protein structure info
                TC = T(:, {'ReactionName', 'Organism', 'GeneID', 'ProteinID', ...
                    'EC Number', 'MetaNetXID','Substrate', 'SMILES', 'InChIKey', 'sequence', 'pdbpath'});
                outFile = fullfile(reconstructionDir, 'kcatData', 'CatPred_input.csv');
                writetable(TC, outFile);
                fprintf('CatPred input written to: %s\n', outFile);
        end
    end

    writtenTable = T;
end

% ------------------------------ Helpers ----------------------------------
function val = safeGetCell(cellarr, idx, defaultVal)
% SAFEGETCELL  Return cell array element if exists and non-empty; else default.
    val = defaultVal;
    try
        if idx >= 1 && idx <= numel(cellarr)
            v = cellarr{idx};
            if ~isempty(v)
                val = v;
            end
        end
    catch
        % fall through with default
    end
end
