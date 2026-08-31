function kcatList_complete = completeKcatMatch(model, OUT, DLmodel, parameters)
% completeKcatMatch (ECOMAP)
%   Build a normalized kcatList struct (rxn, gene, substrate, kcat) from the
%   "complete/common" matches already computed in OUT, mapping ProteinID → gene
%   using the ECOMAP enzymeConstraints layer.
%
%   This utility is robust to slight column-name differences across
%   DL models (DLKcat / UniKP / CatPred) and tables. It:
%     1) Selects the proper OUT.overallTables.* subtable from DLmodel.
%     2) Collects ReactionName, ProteinID, Substrate (with fallbacks), and the
%        target kcat column (prefers experimental kcat if present; otherwise
%        uses the model-specific predicted kcat).
%     3) Maps ProteinID to enzyme index, then to gene IDs via
%        model.enzymeConstraints.{enzymes,genes}.
%     4) Filters out rows with missing rxn IDs or non-positive/non-finite kcat.
%     5) Returns a GECKO/ECOMAP-compatible kcatList structure.
%
% Inputs
%   model     : ECOMAP model, requires model.enzymeConstraints.enzymes & .genes
%   OUT       : struct produced by your matching pipeline; must contain
%               OUT.overallTables.<DLmodelField> as a table.
%               If empty, function will try to load:
%               fullfile(parameters.outputDir, 'AnalyzeKcatMatches.mat')
%               and use OUT_saved.overallTables.
%   DLmodel   : char/string in {'DLKcat','UniKP','CatPred'}
%   parameters: struct with fields:
%               - reconstructionDir (required) : folder holding CSVs (kept for sanity check)
%               - outputDir (optional): where AnalyzeKcatMatches.mat may reside
%
% Output
%   kcatList_complete : struct with fields
%       .source     = 'CompleteMatch'
%       .rxns       : cellstr (ReactionName)
%       .genes      : cellstr (mapped from ProteinID via enzymeConstraints)
%       .substrates : cellstr (preferred original substrate name column)
%       .kcats      : double  (kcat in s^-1; experimental if available, else predicted)
%
% Notes
%   - The function assumes OUT.overallTables.<DLmodel> holds one row per
%     matched (rxn, ProteinID, substrate) entry and contains at least the
%     columns ReactionName and ProteinID, plus substrate and kcat columns
%     under any of the supported aliases below.

    % ----------------------- Parameter checks -----------------------
    if nargin < 4 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.');
        end
    end
    if ~isfield(parameters,'reconstructionDir') || isempty(parameters.reconstructionDir)
        error('parameters.reconstructionDir is required and must point to the folder of CSVs.');
    end
    if ~isfield(parameters,'outputDir') || isempty(parameters.outputDir)
        parameters.outputDir = pwd;
    end

    if nargin < 3 || isempty(DLmodel)
        error('You must specify DLmodel as one of {''DLKcat'',''UniKP'',''CatPred''}.');
    end
    if isstring(DLmodel), DLmodel = char(DLmodel); end
    validModels = {'DLKcat','UniKP','CatPred'};
    if ~ismember(DLmodel, validModels)
        error('DLmodel must be one of {''DLKcat'',''UniKP'',''CatPred''}.');
    end

    % ----------------------- Load OUT if empty ----------------------
    if nargin < 2 || isempty(OUT)
        outfile = fullfile(parameters.outputDir, 'AnalyzeKcatMatches.mat');
        if ~isfile(outfile)
            error('Run the AnalyzeKcatMatches pipeline and save results to AnalyzeKcatMatches.mat.');
        end
        S = load(outfile);
        if ~isfield(S, 'OUT_saved')
            error('AnalyzeKcatMatches.mat does not contain OUT_saved.');
        end
        OUT = S.OUT_saved;
    end

    if ~isfield(OUT, 'overallTables')
        error('OUT.overallTables not found.');
    end

    % Map DLmodel → OUT.overallTables field (tolerate minor naming variations)
    modelField = mapDLField(OUT.overallTables, DLmodel);
    T = OUT.overallTables.(modelField);
    if ~istable(T)
        error('OUT.overallTables.%s is not a table.', modelField);
    end

    % ------------------- Extract required columns -------------------
    rxns       = getCol(T, {'ReactionName','rxn','rxnID','rxns'});
    proteinID  = getCol(T, {'ProteinID','uniprot','Protein_ID','EnzymeID'});
    % Prefer "Substrate" as displayed name; try common fallbacks:
    subs       = getCol(T, {'Substrate','Substrate_name','SubstrateName','Ligand','Reactant','Substrate_norm'}, true);
    % Choose kcat column: prefer experimental if available, else prediction
    if hasAnyCol(T, {'exp_kcat','kcat_exp','experimental_kcat','kcat'})
        kcats = getNumCol(T, {'exp_kcat','kcat_exp','experimental_kcat','kcat'});
    else
        % fallback by DLmodel (column aliases)
        switch DLmodel
            case 'CatPred'
                kcats = getNumCol(T, {'Prediction_(s^(-1))','predicted_kcat','pred_kcat','kcat_pred'});
            otherwise
                kcats = getNumCol(T, {'predicted_kcat','kcat_pred','Prediction_(s^(-1))'});
        end
    end

    % -------------------- Map ProteinID → gene ----------------------
    ec = model.enzymeConstraints;
    req = {'enzymes','genes'};
    for k = 1:numel(req)
        if ~isfield(ec, req{k})
            error('model.enzymeConstraints.%s is missing.', req{k});
        end
    end

    [tf, enzIdx] = ismember(proteinID, ec.enzymes);
    genes = repmat({''}, numel(proteinID), 1);
    genes(tf) = ec.genes(enzIdx(tf));

    % ------------------ Basic quality filtering ---------------------
    ok = ~cellfun(@isempty, rxns) & isfinite(kcats) & (kcats > 0);
    rxns  = rxns(ok);
    genes = genes(ok);
    subs  = subs(ok);
    kcats = kcats(ok);

    % Optional: ensure these reactions exist in enzymeConstraints
    if isfield(ec,'rxns')
        [tf_rxn, ~] = ismember(rxns, ec.rxns);
        if ~all(tf_rxn)
            % keep only those that exist to avoid downstream errors
            rxns  = rxns(tf_rxn);
            genes = genes(tf_rxn);
            subs  = subs(tf_rxn);
            kcats = kcats(tf_rxn);
        end
    end

    % --------------------- Build output struct ----------------------
    kcatList_complete = struct();
    kcatList_complete.source     = 'CompleteMatch';
    kcatList_complete.rxns       = rxns;
    kcatList_complete.genes      = genes;
    kcatList_complete.substrates = subs;
    kcatList_complete.kcats      = kcats;
end

% ====================== helpers (local) ======================

function name = mapDLField(overallTables, DLmodel)
% Tolerant field resolution for OUT.overallTables
    candidates = {DLmodel, lower(DLmodel), upper(DLmodel)};
    fn = fieldnames(overallTables);
    % exact match first
    for i = 1:numel(candidates)
        if isfield(overallTables, candidates{i})
            name = candidates{i}; return;
        end
    end
    % case-insensitive contains
    for i = 1:numel(fn)
        if contains(lower(fn{i}), lower(DLmodel))
            name = fn{i}; return;
        end
    end
    error('Could not resolve OUT.overallTables field for DLmodel=%s', DLmodel);
end

function tf = hasAnyCol(T, names)
    tf = any(ismember(names, T.Properties.VariableNames));
end

function col = getCol(T, names, allowEmpty)
% Return a cellstr column matching the first existing name in `names`.
% If allowEmpty is true and none found, return cell(size(T,1),1) of ''.
    if nargin < 3, allowEmpty = false; end
    found = names(ismember(names, T.Properties.VariableNames));
    if isempty(found)
        if allowEmpty
            col = repmat({''}, height(T), 1);
            return;
        else
            error('Required column missing. Tried: %s', strjoin(names, ', '));
        end
    end
    raw = T.(found{1});
    if iscellstr(raw)
        col = raw;
    elseif isstring(raw)
        col = cellstr(raw);
    elseif iscell(raw)
        % mixed content -> stringify robustly
        col = cellfun(@tochar, raw, 'UniformOutput', false);
    else
        col = cellstr(string(raw));
    end
end

function v = getNumCol(T, names)
% Return a numeric column (double). Accept numeric or string convertible.
    found = names(ismember(names, T.Properties.VariableNames));
    if isempty(found)
        error('Numeric kcat column missing. Tried: %s', strjoin(names, ', '));
    end
    raw = T.(found{1});
    if isnumeric(raw)
        v = double(raw);
    else
        % convert from string/cellstr to double, NaN-safe
        if iscellstr(raw), raw = string(raw); end
        if iscell(raw),    raw = cellfun(@tochar, raw, 'UniformOutput', false); raw = string(raw); end
        v = str2double(raw);
    end
end

function s = tochar(x)
    if ischar(x), s = x; else, s = char(string(x)); end
end
