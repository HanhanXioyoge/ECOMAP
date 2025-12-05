function ecModel = convertecModel(model, ecModeltype, parameters)
% CONVERTECMODEL  Convert a COBRA model into an enzyme-constraint model (ecModel).
%
% ecModel = convertecModel(model, ecModeltype, parameters)
%
% Inputs:
%   model        : COBRA-like model loaded by loadModel (expects model.type == 'TRADITION').
%   ecModeltype  : 'basic' | 'isozyme' | 'integrated'
%                  - basic     : irreversible split; total protein pool bookkeeping
%                  - isozyme   : + isozyme expansion on network
%                  - integrated: + add enzyme pseudo-mets/usage reactions
%   parameters   : struct with fields sigma, Ptot, f, dataDir, uniprot, etc.
%
% Key behavior requested:
%   - Identify "pure compartment switch of the SAME chemical entity" reactions:
%       * exactly 2 metabolites involved;
%       * stoich are opposite and equal in magnitude;
%       * compartments differ;
%       * met Name/Formula/Charge are identical;
%     These reactions are excluded from enzymeConstraints and are also skipped in irreversible splitting and isozyme expansion.
%   - enzymeConstraints.rxns is built after excluding those reactions.
%   - enzymeConstraints.genes is the set of genes that are involved ONLY in the retained
%     (non-transport) reactions; UniProt mapping and rxnEnzMat are built from this gene set.
%
% NOTE:
%   This function builds the ec "skeleton" (protein pool pseudo-met, usage, etc.).
%   Coupling of -1/kcat to catalytic reactions and binding the Ptot*f*sigma budget
%   to 'prot_pool_exchange' should be added in your downstream steps.

    % -------- Parameters --------
    if nargin < 3 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end

    % -------- Select ecModel type (interactive if empty) --------
    if nargin < 2 || isempty(ecModeltype)
        options = {'basic','isozyme','integrated'};
        [indx, tf] = listdlg('PromptString', {'Please select the ecModeltype:'}, ...
                             'SelectionMode','single', ...
                             'ListString',options, ...
                             'Name','Select ecModeltype', ...
                             'ListSize',[320,160]);
        if ~tf, error('No ecModeltype selected.'); end
        ecModeltype = options{indx};
    end

    % -------- Guards on model --------
    if ~isfield(model,'S') || ~isfield(model,'rxns') || ~isfield(model,'mets')
        error('Input "model" must contain S/rxns/mets fields.');
    end
    if ~isfield(model,'type') || ~strcmpi(model.type,'TRADITION')
        error('The input model must be loaded by loadModel and have model.type == ''TRADITION''.');
    end
    
    if isfield(model, 'information')
        model.information.organism              = parameters.org_name;
        model.information.taxonomicID           = parameters.taxonomicID;
        model.information.uniprot.ID            = parameters.uniprot.ID;
        model.information.uniprot.type          = parameters.uniprot.type;
        model.information.uniprot.geneIDfield   = parameters.uniprot.geneIDfield;
        model.information.uniprot.reviewed      = parameters.uniprot.reviewed;
    end

    % -------- Initialize enzymeConstraints header fields --------
    enzymeConstraints = struct();
    enzymeConstraints.sigma = parameters.sigma;
    enzymeConstraints.Ptot  = parameters.Ptot;
    enzymeConstraints.f     = parameters.f;

    % -------- Gene cleanup (remove genes unused) --------
    model = deleteUnusedGenes(model, 0);

    if length(model.genes) ~= length(model.geneMiriams)
        error('The input model needs to be detected and contains insufficient basic information.')
    end

    pseudoGenes = {'s0001'};
    if ismember(pseudoGenes, model.genes)
        model = removemodelGenes(model, pseudoGenes);
    end

    % -------- Step 1: Normalize direction (flip negative-only to positive) -----
    swap_rxns = model.lb < 0 & model.ub == 0;
    if any(swap_rxns)
        model.S(:,swap_rxns) = -model.S(:,swap_rxns);
        model.ub(swap_rxns)  = -model.lb(swap_rxns);
        model.lb(swap_rxns)  = 0;
    end

    % -------- Step 2: Recompute rev flag (LB<0 & UB>0) OR exchange (1 nonzero) --
    model.rev = false(numel(model.rxns),1);
    for i = 1:numel(model.rxns)
        if (model.lb(i) < 0 && model.ub(i) > 0) || nnz(model.S(:,i)) == 1
            model.rev(i) = true;
        end
    end

    % -------- Step 3: Convert to irreversible (exclude exchange + pure transport) ---
    % Pre-detect “pure same-entity compartment switch” reactions so we DO NOT split them.
    [~, exchIdx, ~, pureIdx] = getBoundaryAndTransport(model, 'all');
    
    nRxn = numel(model.rxns);
    toSplitMask = true(nRxn,1);
    toSplitMask(exchIdx) = false;       % exclude exchange/boundary from splitting
    toSplitMask(pureIdx) = false;       % exclude pure-transport from splitting
    
    rxnsToSplit = model.rxns(toSplitMask);
    model = convertToIrrev(model, rxnsToSplit);


    % -------- Step 4: Pre-filtered isozyme expansion (expand ONLY what we will keep) ---
    if strcmp(ecModeltype,'integrated') || strcmp(ecModeltype,'isozyme')
        % Define once and reuse (strict transport filters incl. proton-coupled & stoich-identical)
        opts = struct();
        opts.useNameHeuristic         = true;
        opts.nameHeuristicRequireNoEC = true;
        opts.useProtonCoupling        = true;
        opts.useStoichFallback        = true;
        opts.keywords = { ...
            'transport','transporter','translocase','symport','antiport', ...
            'channel','permease','pump','carrier','porin','shuttle', ...
            'uptake','secretion','efflux','influx','export','import'};
        opts.keepKeywords = { ...
            'abc','pts','phosphotransferase','dehydratase','reductase', ...
            'oxidoreductase','hydrolase','lyase','isomerase','ligase','transferase'};
        opts.extraBlacklist = {};
        opts.extraWhitelist = {};
        opts.protonRegex = '^(h|h\+)\[[a-z]\]$';
        opts.verbose = false;

        % Pass 4a: compute the SAME keep set we will use for enzymeConstraints, on the *current* model
        [pre_keep_idx, pre_report] = buildEnzConstrRxnSet(model, opts);
    
        % Only expand reactions in the pre-keep set
        expandMask = false(numel(model.rxns),1);
        expandMask(pre_keep_idx) = true;
    
        % Expand only the selected reactions (others remain untouched = no wasteful copies)
        [model, ~] = expandModelSelective(model, expandMask);
    
        % Cosmetic: keep reversible/isozyme copies grouped
        model = sortIdentifiers(model);
    else
        model = sortIdentifiers(model);
    end
    
    % -------- Step 5: Build reaction set for enzymeConstraints (AFTER expansion) -------
    if ~isfield(model,'rxnGeneMat')
        error('model.rxnGeneMat is required to build enzymeConstraints.');
    end
    
    % Reuse the exact same opts as Step 4 to stay consistent
    opts = struct();
    opts.useNameHeuristic         = true;
    opts.nameHeuristicRequireNoEC = true;
    opts.useProtonCoupling        = true;
    opts.useStoichFallback        = true;
    opts.keywords = { ...
        'transport','transporter','translocase','symport','antiport', ...
        'channel','permease','pump','carrier','porin','shuttle', ...
        'uptake','secretion','efflux','influx','export','import'};
    opts.keepKeywords = { ...
        'abc','pts','phosphotransferase','dehydratase','reductase', ...
        'oxidoreductase','hydrolase','lyase','isomerase','ligase','transferase'};
    opts.extraBlacklist = {};
    opts.extraWhitelist = {};
    opts.protonRegex = '^(h|h\+)\[[a-z]\]$';
    opts.verbose = false;

    [rxnWG_kept, report] = buildEnzConstrRxnSet(model, opts);
    
    enzymeConstraints.rxns = model.rxns(rxnWG_kept);
    
    % Avoid sparse warnings when printing
    fprintf('enzymeConstraints rxns kept: %d (dropped: %d)\n', ...
        full(double(report.summary.totalKeep)), full(double(report.summary.totalDrop)));


    % -------- Step 6: Determine the GENE set from retained reactions only -------
    % We keep only genes that participate in the retained (non-transport) reactions.
    if isempty(rxnWG_kept)
        error('convertecModel:NoEnzRxns', ...
            ['After exclusions (exchange/boundary & pure-compartment-switch), ', ...
             'no gene-associated reactions remain for enzymeConstraints. ', ...
             'Relax filters or check model annotations/GPRs.']);
    end
    
    geneMaskRetained    = any(model.rxnGeneMat(rxnWG_kept, :), 1);  % 1 x nGenes
    geneIdxKept         = find(geneMaskRetained);
    genesKept           = model.genes(geneIdxKept);
    geneShortNamesKept  = model.geneShortNames(geneIdxKept);
    
    if isempty(geneIdxKept)
        error('convertecModel:NoGenesRetained', ...
            ['No genes are associated with the retained reactions (rxnWG_kept). ', ...
             'Ensure rxnGeneMat/grRules are consistent after expansion and that GPRs exist.']);
    end

    % -------- Step 7: Map retained genes to UniProt (first-match, no scoring) ------
    uniprot_Path = fullfile(parameters.dataDir, 'uniprot.tsv');
    if ~isfile(uniprot_Path)
        DownloadUniProtData(parameters.uniprot, parameters.dataDir);
    end
    dbStruct = ParseUniProtData(uniprot_Path);  % expects .aliasMap, .ID, .MW, .seq, etc.
    
    nG     = numel(genesKept);
    enzIDs = repmat({''}, nG, 1);
    MWs    = nan(nG, 1);
    SEQs   = repmat({''}, nG, 1);
    PDB    = cell(nG,1);
    
    for i = 1:nG
        % Placeholder PDB filename (optional)
        PDB{i,1} = sprintf('seq%d.pdb', i);
    
        % Alias-based exact (case-insensitive) matching:
        % dbStruct.aliasMap maps lowercase alias -> row-index vector in dbStruct
        q = lower(genesKept{i});
        q_2 = lower(geneShortNamesKept{i});

        hits = [];
        if isfield(dbStruct, 'aliasMap') && ~isempty(dbStruct.aliasMap) && isKey(dbStruct.aliasMap, q)
            hits = dbStruct.aliasMap(q);
        end
        
        if isfield(dbStruct, 'aliasMap') && ~isempty(dbStruct.aliasMap) && isKey(dbStruct.aliasMap, q_2)
            hits = dbStruct.aliasMap(q_2);
        end

        % If we have any hit(s), pick the FIRST one (no further selection/scoring)
        if ~isempty(hits)
            pick = hits(1);
            if isfield(dbStruct,'ID')   && numel(dbStruct.ID)   >= pick && ~isempty(dbStruct.ID{pick})
                enzIDs{i} = dbStruct.ID{pick};
            end
            if isfield(dbStruct,'MW')   && numel(dbStruct.MW)   >= pick && ~isnan(dbStruct.MW(pick))
                MWs(i)    = dbStruct.MW(pick);
            end
            if isfield(dbStruct,'seq')  && numel(dbStruct.seq)  >= pick && ~isempty(dbStruct.seq{pick})
                SEQs{i}   = dbStruct.seq{pick};
            end
            continue
        end
    
        % Optional fallback: if aliasMap missing, try a simple contains-split fallback
        if isempty(hits) && isfield(dbStruct,'genes') && ~isempty(dbStruct.genes)
            % Split raw gene field on whitespace and common separators, then compare to q
            % (kept minimal; you already have aliasMap for the robust path)
            for r = 1:numel(dbStruct.genes)
                toks = regexp(lower(dbStruct.genes{r}), '[\s,;/\|]+', 'split');
                toks = toks(~cellfun('isempty', toks));
                if any(strcmp(q, toks))
                    pick = r;
                    if isfield(dbStruct,'ID')   && numel(dbStruct.ID)   >= pick && ~isempty(dbStruct.ID{pick})
                        enzIDs{i} = dbStruct.ID{pick};
                    end
                    if isfield(dbStruct,'MW')   && numel(dbStruct.MW)   >= pick && ~isnan(dbStruct.MW(pick))
                        MWs(i)    = dbStruct.MW(pick);
                    end
                    if isfield(dbStruct,'seq')  && numel(dbStruct.seq)  >= pick && ~isempty(dbStruct.seq{pick})
                        SEQs{i}   = dbStruct.seq{pick};
                    end
                    break
                end
            end
        end
    end
    
    % Fallback: pull UniProt IDs from geneMiriams, ONLY for unmapped genesKept
    unmappedLocal = find(cellfun(@isempty, enzIDs));
    if ~isempty(unmappedLocal) && isfield(model,'geneMiriams')
        needIdxGlobal = geneIdxKept(unmappedLocal);   % indices into model.genes aligned to genesKept(unmappedLocal)
        uniFromM = repmat({''}, numel(needIdxGlobal), 1);
        for k = 1:numel(needIdxGlobal)
            gm = model.geneMiriams{needIdxGlobal(k)};
            if ~isempty(gm)
                idx = find(contains(lower(gm.name), 'uniprot'), 1, 'first');
                if ~isempty(idx) && ~isempty(gm.value{idx})
                    uniFromM{k} = gm.value{idx};
                end
            end
        end
        if any(~cellfun(@isempty, uniFromM))
            [foundIds, MW2, seq2, foundMask] = GetUniProtMWSeqBatch(uniFromM, 'PauseBetweenRequests', 0.2);
            if any(foundMask)
                fillLoc = unmappedLocal(foundMask);   % local positions in genesKept to fill
                for t = 1:numel(fillLoc)
                    loc = fillLoc(t);
                    if ~isempty(foundIds{t}),        enzIDs{loc} = foundIds{t}; end
                    if numel(MW2)  >= t && ~isnan(MW2(t)),     MWs(loc)  = MW2(t);     end
                    if numel(seq2) >= t && ~isempty(seq2{t}),  SEQs{loc} = seq2{t};     end
                end
            end
        end
    end
    
    % Hard guard: if still no UniProt at all, abort early
    if all(cellfun(@isempty, enzIDs))
        error('convertecModel:UniProtMappingEmpty', ...
              'No UniProt IDs could be mapped for retained genes. Check gene names / UniProt DB / alias parsing.');
    end
    
    % Finalize gene-centric fields (ONLY retained genes)
    enzymeConstraints.genes    = genesKept;
    enzymeConstraints.enzymes  = enzIDs;
    enzymeConstraints.mw       = MWs;
    enzymeConstraints.sequence = SEQs;
    enzymeConstraints.PDB      = PDB;

    % -------- Step 8: Build enzymeConstraints.rxns & rxnEnzMat ------------------
    switch ecModeltype
        case {'isozyme','integrated'}
            % Network already expanded (selectively) in Step 4.
            % Keep only retained reactions and build incidence to retained genes.
            enzymeConstraints.rxnIdx = rxnWG_kept;            % rows map to model.rxns(rxnWG_kept)
            enzymeConstraints.rxns   = model.rxns(rxnWG_kept);
    
            nR = numel(rxnWG_kept);
            % Columns initially correspond to ALL retained genes (genesKept)
            % geneIdxKept is the index of retained genes in model.genes
            enzymeConstraints.geneIdx = geneIdxKept(:);        % columns map to model.genes(geneIdxKept)
            nC_all = numel(enzymeConstraints.genes);
    
            % Map global gene index -> local column (0 if not retained)
            geneGlobal2Local = zeros(1, numel(model.genes));
            geneGlobal2Local(geneIdxKept) = 1:nC_all;
    
            % Build sparse rxn x gene incidence (subunit stoichiometry assumed 1)
            rxnEnz = spalloc(nR, nC_all, nnz(model.rxnGeneMat(rxnWG_kept,:)));
            for r = 1:nR
                gIdxGlobal = find(model.rxnGeneMat(rxnWG_kept(r), :) > 0);
                cols = geneGlobal2Local(gIdxGlobal);
                cols = cols(cols > 0);
                if ~isempty(cols)
                    rxnEnz(r, cols) = 1;
                end
            end
    
            % drop gene columns that have NO UniProt mapping (enzymes == '')
            validCols = ~cellfun(@isempty, enzymeConstraints.enzymes);
            if ~all(validCols)
                % Warn (optional)
                % warning('Step8:DropUnmappedGenes','Dropping %d unmapped genes from integrated rxnEnzMat.', sum(~validCols));
                rxnEnz = rxnEnz(:, validCols);
   
                % Slice all gene-centric arrays to keep consistency
                enzymeConstraints.genes    = enzymeConstraints.genes(validCols);
                enzymeConstraints.enzymes  = enzymeConstraints.enzymes(validCols);
                enzymeConstraints.mw       = enzymeConstraints.mw(validCols);
                enzymeConstraints.sequence = enzymeConstraints.sequence(validCols);
                enzymeConstraints.PDB      = enzymeConstraints.PDB(validCols);
                enzymeConstraints.geneIdx  = enzymeConstraints.geneIdx(validCols);
            end
    
            % Store matrix (keep sparse)
            enzymeConstraints.rxnEnzMat = rxnEnz;
    
        case 'basic'
            
            % No network expansion; virtually split isozymes in bookkeeping only.
            % Build prefixed rxn IDs (001_, 002_, ...) for OR-isozymes of retained rxns
            enzymeConstraints.rxns = buildVirtualIsozymeRxnIDs(model, rxnWG_kept);
    
            % Build rxnEnzMat (sum of isozyme copies) with the restricted gene set
            [rxnEnzMat, rowMap] = buildRxnEnzMatBasic(model, rxnWG_kept, enzymeConstraints.genes);
            rxnEnzMat = sparse(rxnEnzMat);   % store as sparse to save memory
            % Row mapping back to model rxn indices (useful for later kcat coupling)
            enzymeConstraints.rxnIdx = rowMap(:);
            
            % Also record geneIdx (columns -> model.genes indices)
            enzymeConstraints.geneIdx = geneIdxKept(:);
            % drop gene columns that have NO UniProt mapping (enzymes == '')
            validCols = ~cellfun(@isempty, enzymeConstraints.enzymes);
            if ~all(validCols)
                % Warn (optional)
                % warning('Step8:DropUnmappedGenes','Dropping %d unmapped genes from integrated rxnEnzMat.', sum(~validCols));
                rxnEnzMat = rxnEnzMat(:, validCols);
   
                % Slice all gene-centric arrays to keep consistency
                enzymeConstraints.genes    = enzymeConstraints.genes(validCols);
                enzymeConstraints.enzymes  = enzymeConstraints.enzymes(validCols);
                enzymeConstraints.mw       = enzymeConstraints.mw(validCols);
                enzymeConstraints.sequence = enzymeConstraints.sequence(validCols);
                enzymeConstraints.PDB      = enzymeConstraints.PDB(validCols);
                enzymeConstraints.geneIdx  = enzymeConstraints.geneIdx(validCols);
            end

            enzymeConstraints.rxnEnzMat = rxnEnzMat;

    end
    % Filter out reactions without enzymes.
    noEnzymeIdx = all(enzymeConstraints.rxnEnzMat == 0, 2);
    enzymeConstraints.rxns(noEnzymeIdx) = [];
    enzymeConstraints.rxnEnzMat(noEnzymeIdx, :) = [];

    % -------- Step 9: Add enzyme pseudo-mets & usage reactions (integrated) -----
    % 9.0 Ensure prot_pool metabolite exists BEFORE creating usage reactions
    poolMetID = 'prot_pool';
    if ~ismember(poolMetID, model.mets)
        pool.mets         = poolMetID;
        pool.metNames     = poolMetID;
        pool.compartments = 'c';
        pool.metNotes     = 'Enzyme-usage protein pool';
        model = addMets(model, pool);
    end
    if strcmp(ecModeltype, 'integrated')
        % Only create pseudo-mets for enzymes with valid UniProt IDs
        validEnzCols = ~cellfun(@isempty, enzymeConstraints.enzymes);
        if any(validEnzCols)
            enzUniprot = enzymeConstraints.enzymes(validEnzCols);
            % uniq list (stable) and per-gene -> prot index mapping
            [uniqUniprot, ~, ic] = unique(enzUniprot, 'stable');  % ic maps each valid column -> uniq index
    
            % 9.1 Add prot_* mets if missing (idempotent)
            protIDsAll   = strcat('prot_', uniqUniprot);
            toAddMask    = ~ismember(protIDsAll, model.mets);
            if any(toAddMask)
                proteinMets.mets         = protIDsAll(toAddMask);
                proteinMets.metNames     = proteinMets.mets;
                proteinMets.compartments = repmat({'c'}, numel(proteinMets.mets), 1);
                if isfield(model,'metMiriams')
                    proteinMets.metMiriams = repmat({struct('name',{{'sbo'}},'value',{{'SBO:0000252'}})}, numel(proteinMets.mets), 1);
                end
                if isfield(model,'metCharges')
                    proteinMets.metCharges = zeros(numel(proteinMets.mets),1);
                end
                proteinMets.metNotes     = repmat({'Enzyme-usage pseudometabolite'}, numel(proteinMets.mets), 1);
                model = addMets(model, proteinMets);
            end
    
            % 9.2 Add usage reactions (idempotent)
            % Convention: {prot_i, prot_pool} = [-1, +1], lb=-1000, ub=0, rev=1
            for k = 1:numel(uniqUniprot)
                protMet = protIDsAll{k};
                usageID = ['usage_' protMet];
                if ~ismember(usageID, model.rxns)
                    usageRxn.rxns         = {usageID};
                    usageRxn.rxnNames     = {usageID};
                    usageRxn.mets         = {{protMet, poolMetID}};
                    usageRxn.stoichCoeffs = {[-1, +1]};
                    usageRxn.lb           = -1000;
                    usageRxn.ub           = 0;
                    usageRxn.rev          = 1;
            
                    % (optional) representative gene for grRules
                    repGene  = '';
                    firstCol = find(validEnzCols);
                    firstCol = firstCol(find(ic == k, 1, 'first'));
                    if ~isempty(firstCol) && firstCol <= numel(enzymeConstraints.genes)
                        repGene = enzymeConstraints.genes{firstCol};
                    end
                    usageRxn.grRules = {''};
            
                    model = addRxns(model, usageRxn);
                end
            end
    
            % 9.3 Save mapping from gene columns -> prot_* indices (for later kcat coupling)
            % Build a vector same length as retained gene columns; 0 means "no UniProt"
            enz2protIdx = zeros(1, numel(enzymeConstraints.genes));
            enz2protIdx(validEnzCols) = ic;  % ic in 1..numel(uniqUniprot) for valid columns
            enzymeConstraints.protMets     = protIDsAll;     % list of unique prot_* ids
            enzymeConstraints.enz2protIdx  = enz2protIdx;    % columns -> prot index (0 if unmapped)
        else
            % No valid UniProt enzymes in integrated mode — abort because we cannot create prot_* usage
            error('convertecModel:NoValidUniProt', ...
                  'Integrated mode requires at least one gene with a valid UniProt ID to create prot_* usage reactions.');
        end
    end
    
    % -------- Step 10: Add protein pool exchange and set budget -------------------
    % Add prot_pool_exchange if missing
    poolExID = 'prot_pool_exchange';
    if ~ismember(poolExID, model.rxns)
        poolRxn.rxns         = poolExID;
        poolRxn.rxnNames     = poolExID;
        poolRxn.mets         = {'prot_pool'};
        poolRxn.stoichCoeffs = {-1};   % negative flux PRODUCES prot_pool up to budget
        % Bounds: allow negative flux by default; will set LB to -(Ptot*f*sigma) below
        poolRxn.lb           = -1000;
        poolRxn.ub           = 0;
        poolRxn.rev          = 1;
        model = addRxns(model, poolRxn);
    end
    
    % Set the protein budget immediately if parameters are available
    if isfield(parameters,'Ptot') && isfield(parameters,'f') && isfield(parameters,'sigma') ...
            && ~isempty(parameters.Ptot) && ~isempty(parameters.f) && ~isempty(parameters.sigma)
        Pbudget = parameters.Ptot * parameters.f * parameters.sigma * 1000;   % units must match your model (e.g., g/gDW)
        % With stoich -1 for prot_pool, LB = -Pbudget allows allocating up to +Pbudget to enzymes via negative usage flux
        model = changeRxnBounds(model, poolExID, -abs(Pbudget), 'l');  % set LB
        model = changeRxnBounds(model, poolExID, 0,             'u');  % keep UB at 0
    else
        % Optional warning if budget not set
        % warning('convertecModel:NoProteinBudget', 'Protein budget (Ptot*f*sigma) not set; prot_pool_exchange left at default bounds.');
    end
    
    % -------- Finalize ----------------------------------------------------------
    fields = fieldnames(enzymeConstraints);
    for i = 1:length(fields)
        if isfield(model.enzymeConstraints, fields{i})
            model.enzymeConstraints.(fields{i}) = enzymeConstraints.(fields{i});
        end
    end
    ecModel         = model;
    ecModel.id      = ['ec',model.id];
    ecModel.type    = ['ECOMAP-',ecModeltype];
    ecModel.enzymeConstraints.ecModeltype = ecModeltype;
end

% =========================================================================
% Helper: build virtual isozyme rxn IDs (for 'basic' mode only)
% Adds 3-digit prefixes 001_, 002_, ... per OR-isozyme unit in grRules.
% =========================================================================
function newRxnIDs = buildVirtualIsozymeRxnIDs(model, rxnIdx)
    orSep = ' or ';
    grRules = model.grRules(rxnIdx);
    orCounts = count(grRules, orSep);
    copies   = orCounts + 1;
    totalCopies = sum(copies);

    newRxnIDs = cell(totalCopies,1);
    rowStart  = cumsum([1; copies(1:end-1)]);
    for k = 1:numel(rxnIdx)
        rID = model.rxns{rxnIdx(k)};
        for j = 1:copies(k)
            newRxnIDs{rowStart(k) + (j-1)} = sprintf('%03d_%s', j, rID);
        end
    end
end

% =========================================================================
% Helper: build rxnEnzMat for 'basic' mode by parsing grRules with OR/AND,
% restricted to a provided gene list (genesKept).
% Returns:
%   rxnEnzMat (sum(copies) x numel(genesKept)), and row mapping if needed.
% =========================================================================
function [rxnEnzMat, rowMap] = buildRxnEnzMatBasic(model, rxnIdx, genesKept)
    orSep  = ' or ';
    grRs   = model.grRules(rxnIdx);
    orCnt  = count(grRs, orSep);
    copies = orCnt + 1;

    nCols  = numel(genesKept);
    nRows  = sum(copies);
    rxnEnzMat = zeros(nRows, nCols);
    rowMap    = zeros(nRows,1);

    rowPtr = 0;
    for k = 1:numel(rxnIdx)
        gpr = grRs{k};
        if isempty(gpr), continue; end
        gpr = strrep(gpr,'(',''); gpr = strrep(gpr,')','');
        gpr = strrep(gpr, orSep, ';');
        isoUnits = regexp(gpr, ';', 'split');   % OR-separated isozyme units
        if isempty(isoUnits), isoUnits = {gpr}; end
        for j = 1:numel(isoUnits)
            rowPtr = rowPtr + 1;
            rowMap(rowPtr) = rxnIdx(k);
            andGenes = regexp(strtrim(isoUnits{j}), ' and ', 'split');
            andGenes = strtrim(andGenes);
            if ~isempty(andGenes)
                loc = ismember(genesKept, andGenes);
                if any(loc)
                    rxnEnzMat(rowPtr, loc) = 1; % subunit stoichiometry = 1
                end
            end
        end
    end
end
