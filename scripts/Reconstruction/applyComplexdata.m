function [model, foundComplex, proposedComplex, appliedComplex] = applyComplexdata(model, complexInfo, UseproposedComplex, parameters)
% APPLYCOMPLEXDATA (match-first, enrich-on-selected-proposal)
% Pass 1: compute exact matches and proposals WITHOUT changing the EC layer.
% Pass 2: apply user-selected proposals (or auto-selection) to EC.rxnEnzMat and model GPR.
%
% UseproposedComplex options:
%   - false (default): do not apply any proposal
%   - true           : auto-apply at most ONE proposal per reaction (original behavior)
%   - 'select'/'ui'  : open an interactive table; apply only checked rows
%   - numeric vector : row indices of proposedComplex to apply
%   - logical vector : length == height(proposedComplex); true rows will be applied
%
% UI table columns (display names):
%   Apply | rxn in model | rxnName | complexID | rxn in complexdata | match (%) | genes | protID_model | protID_complex | stochiometry
%
% Notes:
% - Internally we keep the variable name 'name' (from complexData) for compatibility.
% - When applying, we clear the entire EC.rxnEnzMat row first to avoid leftover stoichiometries.

    % ------------------- Parameters & inputs -------------------
    if nargin < 4 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end
    if nargin < 3 || isempty(UseproposedComplex)
        UseproposedComplex = false;
    end
    if nargin < 2 || isempty(complexInfo)
        complexInfo = getComplexdata(parameters.taxonomicID);
    end

    % decode JSON if needed
    if ischar(complexInfo) || isstring(complexInfo)
        jsonStr     = fileread(complexInfo);
        complexData = jsondecode(jsonStr);
    else
        complexData = complexInfo;
    end

    % guards
    if ~isfield(model,'enzymeConstraints') || isempty(model.enzymeConstraints)
        error('The model does not contain the enzymeConstraints structure.');
    end
    EC = model.enzymeConstraints;
    req = {'rxns','enzymes','rxnEnzMat','ecModeltype'};
    for f = req
        if ~isfield(EC,f{1})
            error('enzymeConstraints.%s is missing.', f{1});
        end
    end



    % ------------------- Load UniProt DB (for gene lookup) -------------------
    if ~isfield(parameters,'reconstructionDir') || isempty(parameters.reconstructionDir)
        error('parameters.reconstructionDir is required to locate uniprot.tsv.');
    end
    uniprot_Path = fullfile(parameters.reconstructionDir, 'uniprot.tsv');
    if ~exist(uniprot_Path,'file')
        error('UniProt TSV not found at: %s', uniprot_Path);
    end
    dbStruct = ParseUniProtData(uniprot_Path);  % expects: ID, genes, geneAliases, MW, seq, ...

    % ------------------- Resolve reaction names for display -------------------
    ecModeltype = EC.ecModeltype;
    switch lower(ecModeltype)
        case {'integrated','isozyme'}
            rxnNamesForDisplay = EC.rxns;
        otherwise
            rxnNamesForDisplay = extractAfter(EC.rxns, 4);   % compatibility with older label style
    end
    nECrxn    = numel(EC.rxns);
    enzNames  = EC.enzymes(:);
    rxnEnzMat = EC.rxnEnzMat;

    % ------------------- Normalize complexData to matrix -------------------
    [protCell, stoichCell] = normalize_complex_entries(complexData);
    nCplx = numel(protCell);

    % Build complex matrix (#complex x #unique_proteins)
    allProt = vertcat(protCell{:});          % long column
    [complexProts, ~, allProtCol] = unique(allProt, 'stable');
    countPerCplx = cellfun(@numel, protCell);
    rows = repelem((1:nCplx)', countPerCplx);
    cols = allProtCol;
    vals = cell2mat(cellfun(@(v) v(:), stoichCell, 'UniformOutput', false));
    complexMatrix = sparse(rows, cols, vals, nCplx, numel(complexProts));

    % Pre-map EC enzymes -> complex matrix columns
    [inComplexProt, enz2col] = ismember(enzNames, complexProts);

    % ------------------- PASS 1: compute exact & proposals (no enrichment) -------------------
    foundRows        = {};
    proposedRows     = {};
    selectedProposals_auto = [];   % rows: [ecRow, complexIdx], one per reaction (auto)
    selectedMask     = false(nECrxn,1);

    for i = 1:nECrxn
        if i > size(rxnEnzMat,1), break; end
        enzMask = rxnEnzMat(i,:) ~= 0;          % enzymes used by this EC reaction
        if ~any(enzMask), continue; end

        modelProts = enzNames(enzMask);
        colIdx     = enz2col(enzMask);
        validMask  = inComplexProt(enzMask) & (colIdx > 0);
        if ~any(validMask), continue; end

        colIdx_valid     = colIdx(validMask);
        modelProts_valid = modelProts(validMask);

        cand = any(complexMatrix(:, colIdx_valid) > 0, 2);
        if ~any(cand), continue; end

        potComplex = find(cand);

        subMatch   = sum(complexMatrix(potComplex, colIdx_valid) > 0, 2);
        totalUnits = sum(complexMatrix(potComplex, :) > 0, 2);
        modUnits   = numel(modelProts);   % denominator uses total enzymes in reaction

        percMatch  = subMatch ./ modUnits;       % coverage (w.r.t. model enzyme count)
        relSize    = totalUnits ./ modUnits;     % complex size relative to model

        % ---- Exact match ----
        exactMask = (percMatch == 1) & (relSize == 1);
        if any(exactMask)
            ii = potComplex(find(exactMask, 1, 'first'));
            [~, modelProtsIdx] = ismember(modelProts_valid, enzNames);
            valsAssign = full(complexMatrix(ii, colIdx_valid));
            rxnEnzMat(i, modelProtsIdx) = valsAssign;

            rxnID      = EC.rxns{i};
            rxnNameVal = get_model_rxnName(model, rxnID);

            foundRows(end+1, :) = { ...
                rxnID, ...
                rxnNameVal, ...
                complexData(ii).complexID, ...
                complexData(ii).name, ... 
                complexData(ii).geneName, ...
                modelProts, ...
                complexData(ii).protID, ...
                complexData(ii).stochiometry ...
            }; %#ok<AGROW>
            continue
        end

        % ---- Proposals (multi-protein only) ----
        if modUnits > 1
            rxnID      = EC.rxns{i};
            rxnNameVal = get_model_rxnName(model, rxnID);

            % A) fully cover model proteins but complex is larger → pick minimal extra size
            moreMask = (percMatch == 1) & (relSize > 1);
            if any(moreMask)
                idxA = find(moreMask);
                [~, kA] = min(relSize(moreMask));
                pickA = potComplex(idxA(kA));
                proposedRows(end+1, :) = { ...
                    rxnID, ...
                    rxnNameVal, ...
                    complexData(pickA).complexID, ...
                    complexData(pickA).name, ...         
                    complexData(pickA).geneName, ...
                    modelProts, ...
                    complexData(pickA).protID, ...
                    complexData(pickA).stochiometry, ...
                    relSize(idxA(kA)) * 100, ...
                    0 ...
                }; %#ok<AGROW>
            end

            % B) partial coverage >= 75% and complex not larger → pick highest coverage
            partMask = (percMatch >= 0.75) & (percMatch < 1) & (relSize <= 1);
            if any(partMask)
                idxB = find(partMask);
                [~, kB] = max(percMatch(partMask));
                pickB = potComplex(idxB(kB));
                proposedRows(end+1, :) = { ...
                    rxnID, ...
                    rxnNameVal, ...
                    complexData(pickB).complexID, ...
                    complexData(pickB).name, ...   
                    complexData(pickB).geneName, ...
                    modelProts, ...
                    complexData(pickB).protID, ...
                    complexData(pickB).stochiometry, ...
                    percMatch(idxB(kB)) * 100, ...
                    0 ...
                }; %#ok<AGROW>
            end

            % Auto-selection (A preferred over B)
            if exist('pickA','var') && ~isempty(pickA)
                selectedProposals_auto(end+1,:) = [i, pickA]; %#ok<AGROW>
                selectedMask(i) = true;
            elseif exist('pickB','var') && ~isempty(pickB)
                selectedProposals_auto(end+1,:) = [i, pickB]; %#ok<AGROW>
                selectedMask(i) = true;
            end
            clear pickA pickB;
        end
    end

    % write back exact matches now
    model.enzymeConstraints.rxnEnzMat = rxnEnzMat;

    ecType = lower(string(model.enzymeConstraints.ecModeltype));
    ec_rxn = model.enzymeConstraints.rxns(:);
    if ~isempty(proposedRows)
        prop_rxn_raw = proposedRows(:,1);
        switch ecType
            case "basic"
                base_rxn     = stripBasicPrefix(ec_rxn);
                proposed_rxn = stripBasicPrefix(prop_rxn_raw);
            otherwise
                base_rxn     = stripIntegratedSuffix(ec_rxn);
                proposed_rxn = stripIntegratedSuffix(prop_rxn_raw);
        end
    
        nProp    = numel(proposed_rxn);
        keepMask = true(nProp,1);
    
        EC   = model.enzymeConstraints;
        Senz = EC.rxnEnzMat;
        Eenz = EC.enzymes(:);
    
        for k = 1:nProp
            key = proposed_rxn{k};
            idx_all_isozyme = find(strcmp(base_rxn, key));
            if isempty(idx_all_isozyme)
                continue;
            end
            protID_complex = proposedRows{k,7};
            protID_complex = asCellStr(protID_complex);  
            for ii = 1:numel(idx_all_isozyme)
                rEC = idx_all_isozyme(ii);
                enzMask = Senz(rEC,:) ~= 0;
                if ~any(enzMask)
                    continue;
                end
                enzymes_used = Eenz(enzMask);
                enzymes_used = asCellStr(enzymes_used);
        
                if numel(enzymes_used) == numel(protID_complex) && ...
                   all(ismember(protID_complex, enzymes_used)) && ...
                   all(ismember(enzymes_used, protID_complex))
        
                    keepMask(k) = false;
                    break;
                end
            end
        end
        proposedRows = proposedRows(keepMask, :);
    end

    % ------------------- Output tables (include rxnName) -------------------
    baseCols = {'rxn','rxnName','complexID','name','genes','protID_model','protID_complex','stochiometry'};
    if isempty(foundRows)
        foundComplex = cell2table(cell(0, numel(baseCols)), 'VariableNames', baseCols);
    else
        foundComplex = cell2table(foundRows, 'VariableNames', baseCols);
    end

    propCols = [baseCols, {'match','applied'}];
    if isempty(proposedRows)
        proposedComplex = cell2table(cell(0, numel(propCols)), 'VariableNames', propCols);
    else
        proposedComplex = cell2table(proposedRows, 'VariableNames', propCols);
        if ~isnumeric(proposedComplex.match);   proposedComplex.match   = double(proposedComplex.match); end
        if ~islogical(proposedComplex.applied); proposedComplex.applied = logical(proposedComplex.applied); end
    end

    % --- NEW (minimal): placeholder for applied rows, same schema as proposedComplex
    appliedComplex = proposedComplex([],:);

    % ------------------- Decide selection mode -------------------
    selection_mode = 'auto';  % 'auto' | 'none' | 'ui' | 'indices'
    indices_to_apply = [];    % rows of proposedComplex

    if islogical(UseproposedComplex) && isscalar(UseproposedComplex)
        if ~UseproposedComplex, selection_mode = 'none'; else, selection_mode = 'auto'; end
    elseif ischar(UseproposedComplex) || isstring(UseproposedComplex)
        s = lower(string(UseproposedComplex));
        if s == "select" || s == "ui" || s == "manual"
            selection_mode = 'ui';
        else
            warning('Unrecognized string for UseproposedComplex: "%s". Fallback to auto.', s);
            selection_mode = 'auto';
        end
    elseif isnumeric(UseproposedComplex)
        selection_mode   = 'indices';
        indices_to_apply = UseproposedComplex(:);
    elseif islogical(UseproposedComplex) && ~isscalar(UseproposedComplex)
        if height(proposedComplex) == numel(UseproposedComplex)
            selection_mode   = 'indices';
            indices_to_apply = find(UseproposedComplex(:));
        else
            warning('Logical selection length mismatch; fallback to auto.');
            selection_mode = 'auto';
        end
    else
        selection_mode = 'auto';
    end

    % ------------------- Build selectedProposals according to mode -------------------
    selectedProposals = [];  % [ecRow, complexIdx]

    if strcmp(selection_mode, 'auto')
        selectedProposals = selectedProposals_auto;

    elseif strcmp(selection_mode, 'none')
        % do nothing

    elseif strcmp(selection_mode, 'ui')
        if height(proposedComplex) == 0
            disp('No proposals to select.');
        else
            mask = prompt_select_proposals_ui(proposedComplex);
            if isempty(mask)
                disp('Selection canceled. Skip applying proposals.');
            else
                indices_to_apply = find(mask);
            end
        end
    end

    % If we have indices_to_apply, convert them to [ecRow, complexIdx]
    if ~isempty(indices_to_apply)
        indices_to_apply = indices_to_apply(indices_to_apply>=1 & indices_to_apply<=height(proposedComplex));
        if ~isempty(indices_to_apply)
            % enforce at most one per rxn (keep highest match)
            indices_to_apply = dedup_keep_best_per_rxn(indices_to_apply, proposedComplex);

            for k = indices_to_apply(:).'
                rxnID  = proposedComplex.rxn{k};
                cplxID = string(proposedComplex.complexID{k});
                % map rxnID -> ecRow
                ecRow = find(strcmp(EC.rxns, rxnID), 1, 'first');
                if isempty(ecRow)
                    rxnID2 = regexprep(rxnID, '^(EC__|EC_|ec__)\s*', '');
                    ecRow  = find(strcmp(EC.rxns, rxnID2), 1, 'first');
                end
                % map complexID -> complex index in complexData
                cIdx  = find(strcmp(string({complexData.complexID}), cplxID), 1, 'first');

                if isempty(ecRow) || isempty(cIdx)
                    warning('Selection row %d cannot be mapped (rxn/ecRow: %s / %d, cplxIdx: %d). Skip.', ...
                            k, rxnID, ternary(isempty(ecRow),-1,ecRow), ternary(isempty(cIdx),-1,cIdx));
                    continue;
                end
                selectedProposals(end+1,:) = [ecRow, cIdx]; %#ok<AGROW>
            end
        end
    end

    % ------------------- PASS 2: apply selection -------------------
    if ~isempty(selectedProposals)
        % Gather missing UniProt required by selected proposals
        haveSet = model.enzymeConstraints.enzymes(:);
        needUni = {};
        for t = 1:size(selectedProposals,1)
            cplxIdx = selectedProposals(t,2);
            uniList = asCellStr(complexData(cplxIdx).protID);
            miss = setdiff(uniList, haveSet);
            needUni = [needUni; miss(:)]; %#ok<AGROW>
        end
        needUni = unique(needUni, 'stable');

        % Enrich EC by adding ONLY the missing UniProt (on demand)
        for k = 1:numel(needUni)
            uniID = needUni{k};
            [geneID, geneName, geneShort, MWval, SEQval] = pick_gene_from_uniprot(dbStruct, uniID);
            try
                [model, ~, ~] = addGeneEC(model, geneID, ...
                    'OnDuplicate','skip', ...
                    'AddEnzymeConstraint', true, ...
                    'UniProtID', uniID, ...
                    'MW', MWval, ...
                    'Sequence', SEQval, ...
                    'PDB','', ...
                    'GeneName', geneName, ...
                    'GeneShortName', geneShort);
            catch ME
                warning('applyComplexdata:addGeneECFailed', ...
                    'Failed to add UniProt %s via addGeneEC: %s. Continue.', uniID, ME.message);
            end
        end

        % Refresh EC references
        EC = model.enzymeConstraints;
        enzNames  = EC.enzymes(:);
        rxnEnzMat = EC.rxnEnzMat;

        % Apply stoichiometries for each selected proposal; then UPDATE model GPR.
        selKey = containers.Map('KeyType','char','ValueType','logical');

        for t = 1:size(selectedProposals,1)
            ecRow = selectedProposals(t,1);
            cIdx  = selectedProposals(t,2);

            rxnID  = EC.rxns{ecRow};
            uniList= asCellStr(complexData(cIdx).protID);
            sVec   = normalizeStoich(complexData(cIdx).stochiometry, numel(uniList)); % 1 x N

            % Map UniProt to EC columns (subset that exist)
            [ok, col] = ismember(uniList, enzNames);
            col  = col(ok);
            vals = sVec(1, ok);
            
            % --- capture previous gene set for this EC row (for basic clause targeting)
            prevGeneIDs = {};
            prevCols = find(rxnEnzMat(ecRow,:) ~= 0);
            if ~isempty(prevCols) && isfield(EC,'genes')
                tmpGenes = EC.genes(prevCols);
                tmpGenes = tmpGenes(:);
                tmpGenes = tmpGenes(~cellfun(@isempty, tmpGenes));
                prevGeneIDs = unique(tmpGenes, 'stable');
            end
            
            % >>> clear the whole row to avoid leftover stoichiometries <<<
            rxnEnzMat(ecRow, :) = 0;
            if ~isempty(col)
                rxnEnzMat(ecRow, col) = vals;
            end

            % >>> clear the whole row to avoid leftover stoichiometries <<<
            rxnEnzMat(ecRow, :) = 0;
            if ~isempty(col)
                rxnEnzMat(ecRow, col) = vals;
            end

            % ------------------- UPDATE MODEL-LEVEL GPR -------------------
            rModel = locate_model_rxn(model.rxns, rxnID);
            if isempty(rModel)
                warning('applyComplexdata:RxnNotFoundInModel', ...
                        'Model.rxns does not contain "%s"; skip GPR update for this reaction.', rxnID);
            else
                % derive newGenes from EC.genes aligned to mapped UniProt columns
                newGenes = {};
                if isfield(EC,'genes')
                    newGenes = EC.genes(col);
                    newGenes = newGenes(:);
                    newGenes = newGenes(~cellfun(@isempty, newGenes));
                    newGenes = unique(newGenes, 'stable');
                end
            
                % ensure grRules / genes / rxnGeneMat ok
                model = ensure_field(model,'grRules',@() repmat({''}, numel(model.rxns),1));
                [model, ~] = ensure_genes_exist_min(model, newGenes);
                model = ensure_rxnGeneMat_size(model);
            
                ecType = '';
                if isfield(EC,'ecModeltype') && ~isempty(EC.ecModeltype)
                    ecType = lower(string(EC.ecModeltype));
                end
            
                if ecType == "basic"
                    % basic：Replace an OR clause → precisely locate with prevGeneIDs; If it doesn't hit, fall back to the numbered prefix/maximum overlap
                    model = replace_basic_clause_and_rebuild_row_targeted(model, rModel, newGenes, prevGeneIDs, rxnID);
                else
                    % isozyme/integrated：
                    model.rxnGeneMat(rModel, :) = 0;
                    if ~isempty(newGenes) && isfield(model,'genes') && ~isempty(model.genes)
                        [~, geneCols] = ismember(newGenes, model.genes);
                        geneCols = geneCols(geneCols>0);
                        if ~isempty(geneCols)
                            model.rxnGeneMat(rModel, unique(geneCols)) = 1;
                        end
                    end
                    if ~isempty(newGenes)
                        model.grRules{rModel} = strjoin(newGenes, ' and ');
                    else
                        model.grRules{rModel} = '';
                    end
                end
            end

            selKey([rxnID,'|',char(complexData(cIdx).complexID)]) = true;
        end

        % mark 'applied' in proposedComplex
        if ~isempty(proposedComplex)
            for u = 1:height(proposedComplex)
                rxnID  = proposedComplex.rxn{u};
                cplxID = char(proposedComplex.complexID{u});
                key = [rxnID,'|',cplxID];
                if isKey(selKey, key)
                    proposedComplex.applied(u) = true;
                end
            end
        end

        % write back updated rxnEnzMat
        model.enzymeConstraints.rxnEnzMat = rxnEnzMat;
        % If propose is used, use the line below to check for excess unused enzymes
        model = cleanupUnusedEnzymes(model);
    end
    
    appliedFromFound = foundComplex(:, baseCols); 
    appliedFromFound.match   = repmat(100, height(foundComplex), 1);
    appliedFromFound.applied = true(height(foundComplex), 1);
    appliedFromFound = appliedFromFound(:, propCols);

    if height(proposedComplex) > 0
        appliedFromProp = proposedComplex(proposedComplex.applied, :);
    else
        appliedFromProp = proposedComplex([],:);
    end

    appliedComplex = [appliedFromFound; appliedFromProp];

    fprintf('Complex mapping: exact=%d, proposed=%d, applied=%d.\n', ...
        height(foundComplex), height(proposedComplex), sum(proposedComplex.applied));

end

% =========================== Local helpers ===========================

function [protCell, stoichCell] = normalize_complex_entries(complexData)
% Normalize protID to col-cellstr and stochiometry to row numeric aligned to protID
    nCplx = numel(complexData);
    protCell   = cell(nCplx,1);
    stoichCell = cell(nCplx,1);

    for i = 1:nCplx
        % protID
        p = complexData(i).protID;
        if iscell(p)
            if any(cellfun(@iscell, p)), p = horzcat(p{:}); end
            p = cellfun(@char, p, 'UniformOutput', false);
        elseif isstring(p) || ischar(p)
            p = {char(p)};
        else
            p = {char(string(p))};
        end
        protCell{i} = p(:);

        % stochiometry
        s = complexData(i).stochiometry;
        if iscell(s)
            if all(cellfun(@isnumeric, s))
                s = cell2mat(s(:))';
            else
                s = strjoin(cellfun(@char, s, 'UniformOutput', false), ',');
            end
        end
        if isstring(s), s = char(s); end
        if ischar(s)
            s2 = strrep(s, ';', ',');
            s2 = strtrim(regexprep(s2, '^\[|\]$', ''));
            toks = regexp(s2, '[, \t]+', 'split');
            nums = str2double(toks);
            if ~isempty(nums) && all(~isnan(nums))
                s = nums(:)';    % row
            else
                s = zeros(1, numel(protCell{i}));
            end
        end
        if isnumeric(s) && all(s==0), s(:) = 1; end
        if iscolumn(s), s = s'; end

        % align length
        if numel(s) < numel(protCell{i})
            s = [s, ones(1, numel(protCell{i}) - numel(s))];
        elseif numel(s) > numel(protCell{i})
            s = s(1:numel(protCell{i}));
        end
        stoichCell{i} = s;
    end
end

function sRow = normalizeStoich(sField, n)
% Convert stochiometry field to row numeric with length n (zeros->ones)
    if iscell(sField)
        if all(cellfun(@isnumeric, sField)), sField = cell2mat(sField(:))'; else, sField = ''; end
    end
    if isstring(sField), sField = char(sField); end
    if ischar(sField)
        s2 = strrep(sField, ';', ',');
        s2 = strtrim(regexprep(s2, '^\[|\]$', ''));
        toks = regexp(s2, '[, \t]+', 'split');
        nums = str2double(toks);
        if ~isempty(nums) && all(~isnan(nums)), sField = nums(:)'; else, sField = zeros(1,n); end
    end
    if isnumeric(sField) && all(sField==0), sField(:) = 1; end
    if iscolumn(sField), sField = sField'; end
    if numel(sField) < n, sField = [sField, ones(1, n-numel(sField))];
    elseif numel(sField) > n, sField = sField(1:n);
    end
    sRow = sField;
end

function C = asCellStr(x)
% Robustly convert to column cellstr
    if iscell(x)
        C = cellfun(@char, x(:), 'UniformOutput', false);
    elseif isstring(x) || ischar(x)
        C = {char(x)};
    else
        C = {char(string(x))};
    end
end

function [geneID, geneName, geneShort, MWval, SEQval] = pick_gene_from_uniprot(dbStruct, uniID)
% Choose a primary gene ID + meta for a UniProt ID
    geneID   = uniID; geneName = ''; geneShort= ''; MWval = NaN; SEQval = '';
    idx = find(strcmp(dbStruct.ID, uniID), 1, 'first');
    if isempty(idx), return; end

    rawGenes = '';
    if isfield(dbStruct,'genes') && numel(dbStruct.genes)>=idx && ~isempty(dbStruct.genes{idx})
        rawGenes = dbStruct.genes{idx};
    end
    primary = '';
    if isfield(dbStruct,'geneAliases') && numel(dbStruct.geneAliases)>=idx && ~isempty(dbStruct.geneAliases{idx})
        al = dbStruct.geneAliases{idx};
        if ~isempty(al), primary = al{1}; end
    end
    if isempty(primary)
        toks = regexp(rawGenes, '\s+', 'split');
        toks = toks(~cellfun(@isempty, toks));
        if ~isempty(toks), primary = toks{1}; end
    end
    if ~isempty(primary)
        geneID   = primary;
        geneShort= primary;
        geneName = rawGenes;
    end
    if isfield(dbStruct,'MW')  && numel(dbStruct.MW)>=idx  && ~isempty(dbStruct.MW(idx)),   MWval  = dbStruct.MW(idx); end
    if isfield(dbStruct,'seq') && numel(dbStruct.seq)>=idx && ~isempty(dbStruct.seq{idx}), SEQval = dbStruct.seq{idx}; end
end

function rModel = locate_model_rxn(modelRxns, rxnID)
% Locate model reaction index by proposed rxn name (try stripping prefixes)
    rModel = find(strcmp(modelRxns, rxnID), 1, 'first');
    if ~isempty(rModel), return; end
    rxnID2 = regexprep(rxnID, '^(?:\d{3}_)+', '');
    rModel = find(strcmp(modelRxns, rxnID2), 1, 'first');
end

function model = ensure_rxnGeneMat_size(model)
% Ensure rxnGeneMat is (#rxns x #genes); grow with zeros if needed.
    if ~isfield(model,'rxnGeneMat') || isempty(model.rxnGeneMat)
        if isfield(model,'genes') && ~isempty(model.genes)
            model.rxnGeneMat = sparse(numel(model.rxns), numel(model.genes));
        else
            model.rxnGeneMat = sparse(numel(model.rxns), 0);
        end
    else
        [rG,cG] = size(model.rxnGeneMat);
        if rG < numel(model.rxns)
            model.rxnGeneMat(numel(model.rxns), max(cG,1)) = 0;   % grow rows
        elseif rG > numel(model.rxns)
            model.rxnGeneMat = model.rxnGeneMat(1:numel(model.rxns), :);
        end
        if isfield(model,'genes') && ~isempty(model.genes)
            if cG < numel(model.genes)
                model.rxnGeneMat(:, numel(model.genes)) = 0;      % grow cols
            elseif cG > numel(model.genes)
                model.rxnGeneMat = model.rxnGeneMat(:, 1:numel(model.genes));
            end
        end
    end
end

function model = ensure_field(model, fld, initFcn)
% Ensure a model-level vector exists with proper length (#rxns)
    if ~isfield(model, fld) || isempty(model.(fld))
        model.(fld) = initFcn();
    else
        if numel(model.(fld)) < numel(model.rxns)
            if iscell(model.(fld))
                model.(fld){numel(model.rxns)} = [];
            else
                model.(fld)(numel(model.rxns),1) = 0;
            end
        elseif numel(model.(fld)) > numel(model.rxns)
            model.(fld) = model.(fld)(1:numel(model.rxns));
        end
    end
end

function mask = prompt_select_proposals_ui(T)
% Show a modal UI table to pick proposals (returns logical mask of size height(T))
% Display columns:
%   Apply | rxn in model | rxnName | complexID | rxn in complexdata | match (%) | genes | protID_model | protID_complex | stochiometry
% Added: a "Select All" checkbox to toggle all Apply states.

    if height(T) == 0
        mask = [];
        return;
    end

    n = height(T);
    Apply = false(n,1);

    % stringify complex fields for display
    genes_str         = cellfun(@(x) join_items(x), T.genes,          'UniformOutput', false);
    prot_model_str    = cellfun(@(x) join_items(x), T.protID_model,   'UniformOutput', false);
    prot_complex_str  = cellfun(@(x) join_items(x), T.protID_complex, 'UniformOutput', false);
    stoich_str        = cellfun(@(x) join_items(x), T.stochiometry,   'UniformOutput', false);

    % compose cell matrix (same row order as T)
    Data = [ ...
        num2cell(Apply), ...
        T.rxn, ...
        T.rxnName, ...
        T.complexID, ...
        T.name, ...          % “rxn in complexdata” (UI header only)
        num2cell(T.match), ...
        genes_str, prot_model_str, prot_complex_str, stoich_str ...
    ];

    colNamesUI = { ...
        'Apply', ...
        'rxn in model', ...
        'rxnName', ...
        'complexID', ...
        'rxn in complexdata', ...
        'match (%)', ...
        'genes', ...
        'protID_model', ...
        'protID_complex', ...
        'stochiometry' ...
    };

    f = uifigure('Name','Select proposals to apply','Position',[100 100 1280 600], 'WindowStyle','modal');
    tbl = uitable(f, ...
        'Data', Data, ...
        'ColumnEditable', [true, false(1, numel(colNamesUI)-1)], ...
        'ColumnName', colNamesUI, ...
        'ColumnWidth', repmat({'auto'},1,numel(colNamesUI)), ...
        'Position', [10 60 1260 500]);

    uilabel(f,'Text','Tick "Apply" for rows to apply. Use "Select All" to toggle all.',...
        'Position',[10 35 700 20]);

    % --- New: "Select All" checkbox (toggles first column for all rows)
    cbAll = uicheckbox(f,'Text','Select All','Position',[270 5 100 25], ...
        'Value', false, ...
        'ValueChangedFcn', @(src,~) onToggleAll(src, tbl));

    % Action buttons
    uibutton(f,'Text','Apply Selected','Position',[10 5 140 25], ...
        'ButtonPushedFcn', @(~,~) uiresume(f));

    uiwait(f);
    if isgraphics(f)
        wasCancel = isappdata(f,'cancel') && getappdata(f,'cancel');
        D = tbl.Data;
        if wasCancel
            mask = [];
        else
            mask = vertcat(D{:,1});   % first column = logical Apply
        end
        delete(f);
    else
        mask = [];
    end

    function onToggleAll(src, tblHandle)
        D = tblHandle.Data;
        if isempty(D), return; end
        fill = logical(src.Value);
        % first column is "Apply"
        D(:,1) = num2cell(repmat(fill, size(D,1), 1));
        tblHandle.Data = D;
    end
end

function s = join_items(x)
% Convert cell/array/char/string to a concise printable string
    if iscell(x)
        y = cellfun(@char, x, 'UniformOutput', false);
        s = strjoin(y,';');
    elseif isnumeric(x)
        if isrow(x); x = x(:); end
        s = sprintf('%g;', x);
        if ~isempty(s); s(end) = []; end
    elseif isstring(x) || ischar(x)
        s = char(x);
    else
        s = char(string(x));
    end
end

function idx_out = dedup_keep_best_per_rxn(indices, T)
% Keep at most one proposal per rxn: choose the one with highest "match".
    if isempty(indices)
        idx_out = indices;
        return;
    end
    rxnList = T.rxn(indices);
    [uRxn,~,g] = unique(rxnList, 'stable');
    idx_out = zeros(numel(uRxn),1);
    for k = 1:numel(uRxn)
        groupIdx = indices(g==k);
        [~,bestPos] = max(T.match(groupIdx));
        idx_out(k) = groupIdx(bestPos);
        if numel(groupIdx) > 1
            warning('Multiple proposals selected for rxn "%s"; keeping the highest match one.', uRxn{k});
        end
    end
end

function y = ternary(cond, a, b)
% tiny helper for inline printing
    if cond, y = a; else, y = b; end
end

function rxnNameVal = get_model_rxnName(model, rxnID)
% Return model.rxnNames{idx} matched by rxnID (with prefix tolerance)
    rxnNameVal = '';
    if isfield(model,'rxnNames') && ~isempty(model.rxnNames)
        rModel = locate_model_rxn(model.rxns, rxnID);
        if ~isempty(rModel) && numel(model.rxnNames) >= rModel && ~isempty(model.rxnNames{rModel})
            rxnNameVal = model.rxnNames{rModel};
        end
    end
end

function [model, geneCols] = ensure_genes_exist_min(model, geneIDs)

    if ~isfield(model,'genes') || isempty(model.genes)
        model.genes = cell(0,1);
    end
    if isempty(geneIDs)
        geneCols = [];
        return;
    end
    [isMem, ~] = ismember(geneIDs, model.genes);
    if any(~isMem)
        model.genes = [model.genes; geneIDs(~isMem)];
        model.rxnGeneMat(:, numel(model.genes)) = 0; % widen
    end
    [~, geneCols] = ismember(geneIDs, model.genes);
end

function model = replace_basic_clause_and_rebuild_row_targeted(model, rModel, newGenes, prevGenes, rxnID)

    gr = model.grRules{rModel};
    clauses = split_or_clauses(gr);
    newSet  = normalize_gene_set(newGenes);
    prevSet = normalize_gene_set(prevGenes);

    idx = pick_clause_index(clauses, prevSet, rxnID);
    if idx < 1 || idx > numel(clauses); idx = 1; end

    clauses{idx} = join_with_and(newSet);

    clauses = dedup_equiv_and_clauses(clauses);

    for i = 1:numel(clauses)
        c = strtrim(clauses{i});
        if contains(c,' and ') && ~(startsWith(c,'(') && endsWith(c,')'))
            clauses{i} = ['(' c ')'];
        else
            clauses{i} = c;
        end
    end
    model.grRules{rModel} = strjoin(clauses, ' or ');

    allGenes = cell(0,1);
    for i = 1:numel(clauses)
        gset = normalize_gene_set(extract_genes_from_clause(clauses{i}));
        allGenes = [allGenes; gset(:)]; %#ok<AGROW>
    end
    allGenes = unique(allGenes, 'stable');
    [model, geneCols] = ensure_genes_exist_min(model, allGenes);
    model = ensure_rxnGeneMat_size(model);

    model.rxnGeneMat(rModel, :) = 0;
    if ~isempty(geneCols)
        model.rxnGeneMat(rModel, unique(geneCols)) = 1;
    end
end

function clauses = split_or_clauses(gr)

    gr = string(gr);
    if strlength(strtrim(gr)) == 0
        clauses = {''};
    else
        c = regexp(char(gr), '\s+or\s+', 'split');
        clauses = c(:)';
    end
end

function idx = pick_clause_index(clauses, prevSet, rxnID)

    idx = -1;
    if ~isempty(prevSet) && ~isempty(clauses)
        best = -inf; arg = 1;
        for i = 1:numel(clauses)
            oldSet = normalize_gene_set(extract_genes_from_clause(clauses{i}));
            ov = numel(intersect(oldSet, prevSet));
            if ov > best
                best = ov; arg = i;
            end
        end
        if best >= 0
            idx = arg;
            return;
        end
    end

    tok = regexp(char(rxnID), '^(?<n>\d{3})_', 'names');
    if ~isempty(tok)
        n = str2double(tok.n);
        if ~isnan(n) && n >= 1 && n <= numel(clauses)
            idx = n; 
            return;
        end
    end
    idx = 1;
end

function clauses_out = dedup_equiv_and_clauses(clauses_in)
    seen = containers.Map('KeyType','char','ValueType','logical');
    out  = {};
    for i = 1:numel(clauses_in)
        S = normalize_gene_set(extract_genes_from_clause(clauses_in{i}));
        key = strjoin(S,'|');  
        if ~isKey(seen, key)
            seen(key) = true;
            out{end+1} = clauses_in{i}; %#ok<AGROW>
        end
    end
    clauses_out = out;
end

function genes = extract_genes_from_clause(clause)
    if isempty(clause)
        genes = {};
        return;
    end
    c = strtrim(regexprep(clause, '^\(|\)$', ''));
    parts = regexp(c, '\s+and\s+', 'split');
    parts = cellfun(@strtrim, parts, 'UniformOutput', false);
    parts = regexprep(parts, '^\(|\)$', '');
    parts = parts(~cellfun(@isempty, parts));
    genes = parts;
end

function out = normalize_gene_set(genes)
    if isempty(genes)
        out = cell(0,1);
        return;
    end
    if ischar(genes) || isstring(genes); genes = cellstr(genes); end
    out = genes(:);
    out = out(~cellfun(@isempty, out));
    out = unique(out, 'stable');
end

function s = join_with_and(genes)
    genes = normalize_gene_set(genes);
    if isempty(genes), s = ''; else, s = strjoin(genes, ' and '); end
end

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
