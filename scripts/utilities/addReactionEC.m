function [model, rxnIdx, geneIdxs, enzIdxs] = addReactionEC(model, rxnID, mets, stoich, lb, ub, varargin)
% addReactionEC
% -------------------------------------------------------------------------
% Add or update a reaction in an enzyme-constrained model (ecModel), and
% optionally register its genes and enzyme usage in the EC layer.
%
% This function:
%   1) Checks that the input model is an enzyme-constrained model
%      (model.enzymeConstraints.ecModeltype).
%   2) Ensures all metabolites in 'mets' exist (creates them via
%      addMetaboliteEC if needed).
%   3) Adds or overwrites a reaction (stoichiometry, bounds, objective,
%      annotations, eccodes) at the main COBRA/RAVEN model level.
%   4) Optionally adds/ensures genes and enzymes via addGeneEC and links
%      the reaction to enzymes in model.enzymeConstraints.rxnEnzMat.
%
% INPUTS:
%   model   : ecModel-like structure with fields:
%             - mets, rxns, S, lb, ub, rev, c, rxnNames, rxnMiriams,
%               eccodes, genes, rxnGeneMat, grRules, enzymeConstraints, ...
%   rxnID   : char/string, reaction ID (stored in model.rxns)
%   mets    : cellstr, list of metabolite IDs (must match model.mets or
%             will be created if missing and MetComps is provided)
%   stoich  : numeric vector, stoichiometric coefficients for 'mets'
%   lb, ub  : numeric scalars, reaction lower and upper bounds
%
% Name-Value pairs (optional):
%   'RxnName'      : char/string, human-readable reaction name
%   'MetComps'     : scalar or list of compartments for missing metabolites
%                    (index or ID; broadcastable to numel(mets))
%   'c'            : numeric scalar, objective coefficient
%   'rxnMiriam'    : struct, MIRIAM annotation for the reaction
%   'ECCodes'      : EC number(s) for the reaction (char/cellstr)
%   'OnDuplicate'  : 'skip' | 'overwrite' | 'error'
%                    policy if rxnID already exists in model.rxns.
%
%   Gene/Enzyme-related (optional):
%   'GeneIDs'        : list of gene IDs (char/string/cellstr)
%   'GeneShortNames' : list of gene short names
%   'GeneNames'      : list of full gene names
%   'UniProtIDs'     : list of UniProt IDs (one per enzyme)
%   'MW'             : list of protein molecular weights (numeric/cell)
%   'Sequence'       : list of protein sequences
%   'PDB'            : list of PDB IDs or strings
%   'ProtComp'       : compartment for protein pseudo-metabolites
%   'GeneRuleMode'   : 'OR' or 'AND' (used if grRules not provided)
%   'grRules'        : custom GPR string; if empty, auto-built from GeneIDs
%
% OUTPUTS:
%   model   : updated ecModel structure
%   rxnIdx  : index of the (new or overwritten) reaction in model.rxns
%   geneIdxs: indices of affected genes in model.genes (if any)
%   enzIdxs : indices of affected enzymes in model.enzymeConstraints.enzymes
%
% -------------------------------------------------------------------------

    % ---------- Guards ----------
    % Ensure the model is an enzyme-constrained model
    assert(isfield(model,'enzymeConstraints') && ...
           isfield(model.enzymeConstraints,'ecModeltype') && ...
           ~isempty(model.enzymeConstraints.ecModeltype), ...
           'addReactionEC:NotEcModel: Missing/empty model.enzymeConstraints.ecModeltype.');

    % Basic argument type checks for mets and stoich
    if ~(iscell(mets) && isnumeric(stoich))
        error('addReactionEC:BadArgs','mets must be cell; stoich numeric.');
    end
    mets   = mets(:);
    stoich = stoich(:);
    if numel(mets) ~= numel(stoich)
        error('addReactionEC:LenMismatch','mets and stoich lengths must match.');
    end

    % Ensure the model has basic metabolite and reaction fields
    if ~isfield(model,'mets') || ~iscell(model.mets)
        error('addReactionEC:BadModel','model.mets must exist.');
    end
    if ~isfield(model,'rxns') || ~iscell(model.rxns)
        model.rxns = cell(0,1);
    end
    if ~isfield(model,'S') || isempty(model.S)
        % Initialize S if missing
        model.S = sparse(numel(model.mets), numel(model.rxns));
    end

    % ---------- Parse options ----------
    % Use inputParser for optional name-value arguments
    p = inputParser; p.FunctionName = 'addReactionEC';
    addRequired(p,'model');
    addRequired(p,'rxnID');
    addRequired(p,'mets');
    addRequired(p,'stoich');
    addRequired(p,'lb');
    addRequired(p,'ub');
    addParameter(p,'RxnName','',@(x)ischar(x)||isstring(x));
    addParameter(p,'MetComps',[],@(x)ischar(x)||isstring(x)||isnumeric(x)||iscell(x));
    addParameter(p,'c',0,@(x)isnumeric(x)&&isscalar(x));
    addParameter(p,'rxnMiriam',[],@(x)isstruct(x)||isempty(x));
    addParameter(p,'ECCodes',[],@(x)ischar(x)||isstring(x)||iscell(x)||isempty(x));
    % NEW: duplicate reaction handling policy
    addParameter(p,'OnDuplicate','skip',@(x)any(strcmpi(x,{'overwrite','skip','error'})));

    % Gene & enzyme options
    addParameter(p,'GeneIDs',[],@(x)ischar(x)||isstring(x)||iscell(x)||isempty(x));
    addParameter(p,'GeneShortNames',[],@(x)ischar(x)||isstring(x)||iscell(x)||isempty(x));
    addParameter(p,'GeneNames',[],@(x)ischar(x)||isstring(x)||iscell(x)||isempty(x));
    addParameter(p,'UniProtIDs',[],@(x)ischar(x)||isstring(x)||iscell(x)||isempty(x));
    addParameter(p,'MW',NaN,@(x)isnumeric(x)||iscell(x));
    addParameter(p,'Sequence','',@(x)ischar(x)||isstring(x)||iscell(x));
    addParameter(p,'PDB','',@(x)ischar(x)||isstring(x)||iscell(x));
    addParameter(p,'ProtComp','c',@(x)ischar(x)||isstring(x)||isnumeric(x));
    addParameter(p,'GeneRuleMode','OR',@(x)ischar(x)||isstring(x));
    addParameter(p,'grRules','',@(x)ischar(x)||isstring(x));
    parse(p, model, rxnID, mets, stoich, lb, ub, varargin{:});

    % Extract parsed options into local variables
    rxnID     = char(p.Results.rxnID);
    rxnName   = char(p.Results.RxnName);
    metComps  = p.Results.MetComps;
    ccoef     = p.Results.c;
    rxnMir    = p.Results.rxnMiriam;
    ECC       = p.Results.ECCodes;
    dupMode   = lower(char(p.Results.OnDuplicate));   % 'overwrite' | 'skip' | 'error'

    % Gene / enzyme-related parsed values
    GeneIDs   = toCell(p.Results.GeneIDs);
    GShorts   = toCell(p.Results.GeneShortNames);
    GNames    = toCell(p.Results.GeneNames);
    UniIDs    = toCell(p.Results.UniProtIDs);
    MWs       = toCellnum(p.Results.MW);
    SEQs      = toCell(p.Results.Sequence);
    PDBs      = toCell(p.Results.PDB);
    ProtComp  = p.Results.ProtComp;
    ruleMode  = upper(char(p.Results.GeneRuleMode));  % 'OR' or 'AND'
    grRuleIn  = char(p.Results.grRules);

    % Determine reaction reversibility from bounds
    rev = (lb < 0) && (ub > 0);

    % ---------- Duplicate policy at MODEL level (prework) ----------
    % If reaction ID already exists, handle according to OnDuplicate
    existIdx = find(strcmp(model.rxns, rxnID), 1, 'first');
    if ~isempty(existIdx)
        switch dupMode
            case 'skip'
                % Return without changing the model
                warning('addReactionEC:DuplicateModelRxn','Reaction "%s" exists; skip unchanged.', rxnID);
                rxnIdx = existIdx; geneIdxs = []; enzIdxs = [];
                return
            case 'error'
                % Hard error
                error('addReactionEC:DuplicateModelRxn','Reaction "%s" already exists.', rxnID);
            case 'overwrite'
                % Fall through: we will overwrite in-place later
        end
    end

    % ---------- Ensure rxn-level vectors aligned ----------
    % Make sure that all core reaction attributes exist and match nRxns
    nRxns = numel(model.rxns);
    nMets = numel(model.mets);
    model = ensureField(model,'rxnNames',  @() repmat({''}, nRxns,1));
    model = ensureField(model,'lb',        @() zeros(nRxns,1));
    model = ensureField(model,'ub',        @() zeros(nRxns,1));
    model = ensureField(model,'rev',       @() false(nRxns,1));
    model = ensureField(model,'c',         @() zeros(nRxns,1));
    model = ensureField(model,'rxnMiriams',@() cell(nRxns,1));

    % eccodes: store EC numbers per reaction
    if ~isfield(model,'eccodes') || isempty(model.eccodes)
        model.eccodes = repmat({''}, nRxns,1);
    elseif numel(model.eccodes) ~= nRxns
        error('addReactionEC:LenMismatch','model.eccodes length mismatch.');
    end

    % Empty default MIRIAM annotation for reactions
    rxnMirEmpty = struct('name',{},'value',{});
    if isempty(rxnMir), rxnMir = rxnMirEmpty; end

    % ---------- 1) Ensure metabolites exist ----------
    % For each metabolite in 'mets', check if it exists in model.mets.
    % If not, and MetComps is provided, call addMetaboliteEC to create it.
    metIdxs = zeros(numel(mets),1);
    missing = false(numel(mets),1);
    for k = 1:numel(mets)
        hit = find(strcmp(model.mets, mets{k}), 1, 'first');
        if isempty(hit)
            missing(k) = true;
        else
            metIdxs(k) = hit;
        end
    end

    if any(missing)
        % For missing metabolites we need compartment information
        if isempty(metComps)
            error('addReactionEC:MetCompRequired', ...
                  'MetComps required for missing metabolites: %s', strjoin(mets(missing),', '));
        end
        % If single comp provided, broadcast to all mets
        if ischar(metComps) || isstring(metComps) || isnumeric(metComps)
            metComps = repmat({metComps}, numel(mets),1);
        else
            metComps = metComps(:);
            if numel(metComps) ~= numel(mets)
                error('addReactionEC:MetCompsLen','MetComps length must match mets (or be scalar/broadcastable).');
            end
        end

        % Ensure addMetaboliteEC is on the path
        assert(exist('addMetaboliteEC','file')==2 || exist('addMetaboliteEC','file')==6, ...
               'addReactionEC:MissingAddMetaboliteEC: addMetaboliteEC.m not found.');

        % Create each missing metabolite and store its index
        for k = 1:numel(mets)
            if missing(k)
                model = addMetaboliteEC(model, mets{k}, metComps{k});
                metIdxs(k) = numel(model.mets); % newly appended index
            end
        end
        nMets = numel(model.mets);
        % After adding metabolites, S must have one row per metabolite
        if size(model.S,1) ~= nMets
            error('addReactionEC:SRows','S rows must equal numel(model.mets) after adding metabolites.');
        end
    end

    % ---------- 2) Ensure genes/enzymes if enzymatic ----------
    % If GeneIDs or UniProtIDs are provided, treat this as an enzymatic
    % reaction and ensure corresponding genes/enzyme entries exist.
    isEnzymatic = ~isempty(GeneIDs) || ~isempty(UniIDs);
    geneIdxs = [];
    enzIdxs  = [];
    if isEnzymatic
        % addGeneEC must be available
        assert(exist('addGeneEC','file')==2 || exist('addGeneEC','file')==6, ...
               'addReactionEC:MissingAddGeneEC: addGeneEC.m not found.');

        % Normalize and broadcast gene/enzyme lists to a consistent length
        [GeneIDs, GShorts, GNames, UniIDs, MWs, SEQs, PDBs] = ...
            alignGeneEnzLists(GeneIDs, GShorts, GNames, UniIDs, MWs, SEQs, PDBs);

        % For each gene (and optionally enzyme), call addGeneEC
        for i = 1:numel(GeneIDs)
            args = {'GeneName',GNames{i}, 'GeneShortName',GShorts{i}, 'EnsureUnique',true};
            if ~isempty(UniIDs{i})
                % If UniProtID is given, also add enzyme constraints
                args = [args, {'AddEnzymeConstraint', true, ...
                               'UniProtID', UniIDs{i}, 'MW', MWs{i}, ...
                               'Sequence', SEQs{i}, 'PDB', PDBs{i}, ...
                               'ProtComp', ProtComp}];
            end
            [model, gIdx, eIdx] = addGeneEC(model, GeneIDs{i}, args{:});
            geneIdxs(end+1) = gIdx; %#ok<AGROW>
            if ~isempty(eIdx), enzIdxs(end+1) = eIdx; end %#ok<AGROW>
        end

        % If we expected enzymes but none were registered, warn/error
        if isempty(enzIdxs)
            error('addReactionEC:NoEnzymeIdx', ...
                  'Enzymatic reaction requested but no enzymes were registered (provide UniProtIDs or pre-map).');
        end
    end

    % ---------- 3) ADD or REPLACE at MODEL level ----------
    if isempty(existIdx)
        % --- Append a new reaction ---
        [rS,cS] = size(model.S);
        if rS ~= nMets
            error('addReactionEC:SDimMismatch','size(S,1) != numel(model.mets).');
        end

        % Extend S with one new (initially zero) column
        model.S = [model.S, sparse(nMets,1)];
        newCol  = cS + 1;

        % Write stoichiometry into the new column
        for k = 1:numel(metIdxs)
            model.S(metIdxs(k), newCol) = stoich(k);
        end

        % Append reaction-level attributes
        model.rxns      = [model.rxns; {rxnID}];
        model.rxnNames  = [model.rxnNames; { ternary(~isempty(rxnName), rxnName, rxnID) }];
        model.lb        = [model.lb; lb];
        model.ub        = [model.ub; ub];
        model.rev       = [model.rev; rev];
        model.c         = [model.c;  ccoef];
        model.rxnMiriams= [model.rxnMiriams; {rxnMir}];
        model.eccodes   = [model.eccodes; { toECC(ECC) }];
        rxnIdx = numel(model.rxns);

        % Ensure rxnGeneMat exists and is aligned, then set gene association
        model = ensureGrGeneMats(model, rxnIdx);
        if ~isempty(geneIdxs)
            model.rxnGeneMat(rxnIdx, unique(geneIdxs)) = 1;
        end

        % Ensure grRules exists and append a new entry (buildRule if needed)
        model = ensureField(model,'grRules', @() repmat({''}, rxnIdx-1,1));
        model.grRules = [model.grRules; { buildRule(model, geneIdxs, grRuleIn, ruleMode) }];

    else
        % --- Overwrite an existing reaction in-place ---
        rxnIdx = existIdx;

        % Reset stoichiometry column and write new coefficients
        model.S(:, rxnIdx) = 0;
        for k = 1:numel(metIdxs)
            model.S(metIdxs(k), rxnIdx) = stoich(k);
        end

        % Replace reaction attributes
        model.rxnNames{rxnIdx}   = ternary(~isempty(rxnName), rxnName, rxnID);
        model.lb(rxnIdx,1)       = lb;
        model.ub(rxnIdx,1)       = ub;
        model.rev(rxnIdx,1)      = rev;
        model.c(rxnIdx,1)        = ccoef;
        model.rxnMiriams{rxnIdx} = rxnMir;
        model.eccodes{rxnIdx}    = toECC(ECC);

        % Ensure gene matrices exist and reset this reaction's row
        model = ensureGrGeneMats(model, numel(model.rxns));
        model.rxnGeneMat(rxnIdx, :) = 0;
        if ~isempty(geneIdxs)
            model.rxnGeneMat(rxnIdx, unique(geneIdxs)) = 1;
        end

        % Ensure grRules exists and update this reaction's GPR rule
        model = ensureField(model,'grRules', @() repmat({''}, numel(model.rxns),1));
        model.grRules{rxnIdx} = buildRule(model, geneIdxs, grRuleIn, ruleMode);
    end

    % ---------- 4) EC layer (ADD/REPLACE/ SKIP) ----------
    % Update enzymeConstraints.rxns and rxnEnzMat so that the new reaction
    % is connected to the enzyme(s) (by enzyme index).
    if isEnzymatic
        EC = model.enzymeConstraints;

        % Initialize EC.enzymes and EC.rxns if they are empty
        if ~isfield(EC,'enzymes') || isempty(EC.enzymes)
            EC.enzymes = cell(0,1);
        end
        if ~isfield(EC,'rxns')    || isempty(EC.rxns)
            EC.rxns    = cell(0,1);
        end
        nEnz = numel(EC.enzymes);

        % Ensure EC.rxnEnzMat exists and has correct size:
        %   rows   = number of EC.rxns
        %   columns= number of EC.enzymes
        if ~isfield(EC,'rxnEnzMat') || isempty(EC.rxnEnzMat)
            EC.rxnEnzMat = sparse(numel(EC.rxns), nEnz);
        else
            [rE,cE] = size(EC.rxnEnzMat);
            % Extend columns if new enzymes were added
            if cE < nEnz
                EC.rxnEnzMat(:, cE+1:nEnz) = 0;
            elseif cE > nEnz
                error('addReactionEC:ECMatCols','EC.rxnEnzMat columns > #enzymes.');
            end
            % Check row alignment with EC.rxns
            if rE ~= numel(EC.rxns)
                error('addReactionEC:ECMatRows','EC.rxnEnzMat rows must equal numel(EC.rxns) BEFORE edit.');
            end
        end

        % Find existing EC row for this reaction, if any
        ecRow = find(strcmp(EC.rxns, rxnID), 1, 'first');
        if ~isempty(ecRow)
            % Reaction already present in EC layer, respect duplicate policy
            switch dupMode
                case 'skip'
                    warning('addReactionEC:DuplicateECRxn','EC.rxns already has "%s"; skip EC mapping unchanged.', rxnID);
                    model.enzymeConstraints = EC;
                    return
                case 'error'
                    error('addReactionEC:DuplicateECRxn','EC.rxns already has "%s".', rxnID);
                case 'overwrite'
                    % Reset row and continue to set new enzyme mapping
                    EC.rxnEnzMat(ecRow, :) = 0;
            end
        else
            % Append new EC row for this reaction
            EC.rxns = [EC.rxns; {rxnID}];
            EC.rxnEnzMat = [EC.rxnEnzMat; sparse(1, nEnz)];
            ecRow = numel(EC.rxns);
        end

        % Set EC.rxnEnzMat row to point to the relevant enzymes
        if ~isempty(enzIdxs)
            EC.rxnEnzMat(ecRow, unique(enzIdxs)) = 1;
        end
        model.enzymeConstraints = EC;
    end
end

% ================= Helpers =================

function model = ensureField(model, fld, initFcn)
% ensureField
% Ensure a field exists and is non-empty; if missing or empty, initialize
% it by calling the provided initFcn().
    if ~isfield(model, fld) || isempty(model.(fld))
        model.(fld) = initFcn();
    end
end

function c = toCell(x)
% toCell
% Normalize char/string/cellstr or empty input into a column cell array.
    if isempty(x), c = {}; return; end
    if iscell(x),   c = x(:); return; end
    if isstring(x), c = cellstr(x(:)); return; end
    if ischar(x),   c = {x}; return; end
    error('addReactionEC:TypeMismatch','Expected char/string/cellstr.');
end

function c = toCellnum(x)
% toCellnum
% Normalize numeric or cell input into a cell array of numeric values.
    if iscell(x), c = x(:); return; end
    if isnumeric(x)
        if isscalar(x)
            c = {x};
        else
            c = num2cell(x(:));
        end
        return;
    end
    if isempty(x), c = {NaN}; return; end
    error('addReactionEC:TypeMismatch','Expected numeric or cell.');
end

function [G,GS,GN,U,M,S,P] = alignGeneEnzLists(G,GS,GN,U,M,S,P)
% alignGeneEnzLists
% Broadcast and align gene/enzyme-related lists so that they all have the
% same length nG (number of genes). If G is empty but UniProt IDs U exist,
% we use U as fallback gene IDs.
    nG = numel(G);
    if nG==0 && numel(U)>0
        G = U;
        nG = numel(U);
    end
    GS = broadcastTo(GS, nG, '');  % gene short names
    GN = broadcastTo(GN, nG, '');  % gene long names
    U  = broadcastTo(U,  nG, '');  % UniProt IDs
    M  = broadcastTo(M,  nG, NaN);% MW
    S  = broadcastTo(S,  nG, '');  % sequences
    P  = broadcastTo(P,  nG, '');  % PDB IDs
end

function out = broadcastTo(v, N, fill)
% broadcastTo
% Utility to broadcast or fill a list v to length N.
% - If v is cell:
%       * empty -> fill with {fill}
%       * scalar -> replicate N times
%       * length N -> return as column
%       * other -> error
% - If v is numeric scalar -> replicate and wrap in cell
% - If v is char/string -> replicate {char(v)} N times
% - If empty -> fill with {fill}
    if iscell(v)
        if isempty(v)
            out = repmat({fill}, N,1);
        elseif numel(v)==1
            out = repmat(v, N,1);
        elseif numel(v)==N
            out = v(:);
        else
            error('addReactionEC:BroadcastLen','List length mismatch.');
        end
        return;
    end
    if isnumeric(v) && isscalar(v)
        out = num2cell(repmat(v,N,1));
        return;
    end
    if ischar(v) || isstring(v)
        out = repmat({char(v)}, N,1);
        return;
    end
    if isempty(v)
        out = repmat({fill}, N,1);
        return;
    end
    error('addReactionEC:BroadcastType','Unsupported broadcast type.');
end

function y = ternary(cond, a, b)
% ternary
% Simple ternary operator: returns a if cond is true, otherwise b.
    if cond, y = a; else, y = b; end
end

function s = toECC(ECC)
% toECC
% Normalize ECCodes input to a single string (semicolon-separated if cell).
    if isempty(ECC), s = ''; return; end
    if iscell(ECC)
        s = strjoin(cellfun(@char,ECC,'UniformOutput',false), ';');
    else
        s = char(ECC);
    end
end

function rule = buildRule(model, geneIdxs, grRuleIn, mode)
% buildRule
% Build a GPR rule string given gene indices and a join mode ('OR'/'AND'),
% unless a custom grRuleIn is provided, in which case that is used.
    if ~isempty(grRuleIn)
        rule = grRuleIn;
        return;
    end
    if isempty(geneIdxs)
        rule = '';
        return;
    end
    gIDs = model.genes(unique(geneIdxs));
    % Join with ' and ' or ' or ' depending on mode
    rule = strjoin(gIDs, ternary(strcmpi(mode,'AND'),' and ',' or '));
end

function model = ensureGrGeneMats(model, nRxnsNow)
% ensureGrGeneMats
% Make sure model.rxnGeneMat exists and is consistent with model.genes and
% the number of reactions nRxnsNow-1 (when appending a new reaction).
% This function only ensures dimensions, and does not fill in gene links.
    nGenes = isfield(model,'genes') * numel(model.genes);

    if ~isfield(model,'rxnGeneMat') || isempty(model.rxnGeneMat)
        model.rxnGeneMat = sparse(nRxnsNow, nGenes);
    end

    [rG,cG] = size(model.rxnGeneMat);

    if cG < nGenes
        model.rxnGeneMat(:, cG+1:nGenes) = 0;
    elseif cG > nGenes
        error('addReactionEC:rxnGeneMatCols','rxnGeneMat has more columns than model.genes.');
    end

    if rG < nRxnsNow
        model.rxnGeneMat(nRxnsNow, 1) = 0;
    elseif rG > nRxnsNow
        error('addReactionEC:rxnGeneMatRows','rxnGeneMat rows inconsistent.');
    end
end
