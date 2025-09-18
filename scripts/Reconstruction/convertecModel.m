function ecModel = convertecModel(model, ecModeltype, parameters)
% CONVERTECMODEL Converts a COBRA model into an enzyme‐constraint model (ecModel).
%
% ecModel = convertEcModel(model, ecModeltype) converts the input COBRA model
% into an ecModel structure based on the desired ecModeltype.
%
% Input:
%     model       - A COBRA model load by function loadModel
%                   which must have model.name = 'COBRA'.
%     ecModeltype - A string specifying the type of ecModel to build. Three options:
%                     'basic'      : Reversible reactions splitting;Total protein constraints
%                     'isozyme'    : Reversible reactions splitting;Total protein constraints;Isozyme reactions splitting
%                     'integrated' : Reversible reactions splitting;Total protein constraints;Isozymes splitting;Enzymes usage adding
% Output:
%     ecModel     - The converted enzyme‐constraint model.
% Both ecModel types include a protein pool reaction, which modifies the stoichiometric
% matrix S. If model.name is not 'COBRA', the function raises an error.


    % Load parameters once
    if nargin < 3 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.')
        end
    end

    % Select the ecModeltype you want to transform
    if nargin < 2 || isempty(ecModeltype)
        options = {'basic', 'isozyme', 'integrated'};
        [indx, tf] = listdlg('PromptString', {'Please select the ecModeltype:'}, ...
                             'SelectionMode', 'single', ...
                             'ListString', options, ...
                             'Name', 'Select ecModeltype', ...
                             'ListSize', [320, 160]);
        if ~tf
            error('No ecModeltype selected. Operation cancelled.');
        end
        ecModeltype = options{indx};
    end

    % Check Input Model Type
    if ~isfield(model, 'name') || ~strcmpi(model.name, 'TRADITION')
        error('The input model must be loaded using the loadModel function.');
    end

    enzymeConstraints.sigma = parameters.sigma; 
    enzymeConstraints.Ptot = parameters.Ptot;
    enzymeConstraints.f = parameters.f;

    % Remove pseudo-genes and unusedGenes
    model = deleteUnusedGenes(model, 0);
    pseudoGenes = {};
    for j = 1:length(model.genes)
        geneMiriam = model.geneMiriams{j};
        if isempty(geneMiriam)
            pseudoGenes{end+1} = model.genes{j};
        end
    end

    if ~isempty(pseudoGenes)
        model = removeGenes(model, pseudoGenes);
    end

    % Step 1: Reverse the direction of reactions that have been defined to have a negative flux only.
    to_swap=model.lb < 0 & model.ub == 0;
    model.S(:,to_swap)=-model.S(:,to_swap);
    model.ub(to_swap)=-model.lb(to_swap);
    model.lb(to_swap)=0;

    % Step 2: Correct rev vector: true if LB < 0 & UB > 0, or it is an exchange reaction:
    model.rev = false(size(model.rxns));
    for i = 1:length(model.rxns)
        if (model.lb(i) < 0 && model.ub(i) > 0) || sum(model.S(:,i) ~= 0) == 1
            model.rev(i) = true;
        end
    end

    % Step 3: Make irreversible model (appends _REV to reaction IDs to indicate reverse reactions)
    [~,exchRxns] = getExchangeRxns(model);
    nonExchRxns = model.rxns;
    nonExchRxns(exchRxns) = [];
    model=convertToIrrev(model, nonExchRxns);

    % Step 4: Expand model, to separate isozymes (appends _EXP_* to reaction IDs to indicate duplication)
    if strcmp(ecModeltype, 'integrated') || strcmp(ecModeltype, 'isozyme')
        model=expandModel(model);
        % Sort reactions, so that reversible and isozymic reactions are kept near.
        model=sortIdentifiers(model);
    end

    % Step 5: Make enzymeConstraints structure, one for gene-associated reaction.
    rxnWithGene  = find(sum(model.rxnGeneMat,2));

    if strcmp(ecModeltype, 'integrated') || strcmp(ecModeltype, 'isozyme')
        enzymeConstraints.rxns = model.rxns(rxnWithGene);    

    elseif strcmp(ecModeltype, 'basic')
        orSeparator   = ' or ';      % Separator used in the GPR string to denote "or"
        prefixFormat  = '%03d_';     % Prefix format for isozyme numbering (e.g., '001_')
        maxIsozymes   = 999;         % Maximum allowed number of isozymes for one reaction

        %Calculate the number of "or" occurrences in the GPRs for reactions with genes,
        %and determine the number of copies (isozymes) for each reaction.
        orCounts = count(model.grRules(rxnWithGene), orSeparator);
        copies   = orCounts + 1;

        %Create an array of indices by repeating each reaction index by its copy number.
        numRxnsOriginal = length(orCounts);
        copiedIndices   = repelem(rxnWithGene, copies).';

        %Retrieve the reaction IDs from the model corresponding to the copied indices.
        tempRxnIDs = model.rxns(copiedIndices);
        newRxnIDs  = tempRxnIDs;

        %Loop over the original reaction list to add a numbered prefix to each isozyme copy.
        nextIndex = 1;
        for i = 1:numel(model.rxns)
            isozymeCounter = 1;
            % Check if the current reaction matches the next reaction in the temporary list.
            if nextIndex <= length(tempRxnIDs) && strcmp(model.rxns{i}, tempRxnIDs{nextIndex})
                while true
                    %Generate prefix string using the specified format.
                    prefixStr = compose(prefixFormat, isozymeCounter);
                    %Concatenate the prefix with the original reaction ID.
                    newRxnIDs{nextIndex} = [prefixStr{1} tempRxnIDs{nextIndex}];

                    isozymeCounter = isozymeCounter + 1;
                    %Throw an error if the maximum allowed isozyme count is exceeded.
                    if isozymeCounter > maxIsozymes
                        error('Increase index size: maximum isozyme count of %d exceeded for reaction %s.', maxIsozymes, model.rxns{i});
                    end

                    nextIndex = nextIndex + 1;
                    %Exit the loop if all copies of the current reaction are processed.
                    if nextIndex > length(tempRxnIDs) || ~strcmp(model.rxns{i}, tempRxnIDs{nextIndex})
                        break;
                    end
                end
            end
        end
        enzymeConstraints.rxns = newRxnIDs;
    end

    % Step 6: Gather enzyme information
    % Download proteomic data from the Uniprot database
    uniprot_Path = fullfile(parameters.dataDir, 'uniprot.tsv');
    if ~isfile(uniprot_Path)
        DownloadUniProtData(parameters.uniprot, parameters.dataDir);
    end
    dbStruct = ParseUniProtData(uniprot_Path);

    Numl    = numel(model.genes);
    Lia     = false(Numl,1);
    Loc     = zeros(Numl,1);
    PDB     = {};

    for i = 1:Numl
        PDB{end+1,1}  = sprintf('seq%d.pdb', i);
        matches = cellfun(@(x) contains(x, model.genes{i}), dbStruct.genes);
        if any(matches)
            Lia(i) = true;
            Loc(i) = find(matches, 1); 
        end
    end
    
    % Build matched enzymeConstraints (only indices where Lia is true)
    matchedIdx                      = find(Lia);
    enzymeConstraints.genes         = model.genes(matchedIdx);    
    enzymeConstraints.enzymes       = dbStruct.ID(Loc(matchedIdx));
    enzymeConstraints.mw            = dbStruct.MW(Loc(matchedIdx));
    enzymeConstraints.sequence      = dbStruct.seq(Loc(matchedIdx)); 
    enzymeConstraints.PDB           = PDB;

    % Build Nomatch structure and try to extract UniProt IDs from model.geneMiriams
    NoMatchGenes = ~Lia;
    if any(NoMatchGenes)
        NoMatchIdx   = find(NoMatchGenes);
        Nomatch.Genes = model.genes(NoMatchGenes);
        Nomatch.enzymes = cell(numel(NoMatchIdx),1); % keep 1:1 correspondence
        
        for k = 1:numel(NoMatchIdx)
            i = NoMatchIdx(k);
            geneMiriam = model.geneMiriams{i};
            if ~isempty(geneMiriam)
                idx = find(contains(lower(geneMiriam.name), 'uniprot'), 1, 'first');
                if ~isempty(idx) && ~isempty(geneMiriam.value{idx})
                    Nomatch.enzymes{k} = geneMiriam.value{idx}; % put the ID string or cell
                else
                    Nomatch.enzymes{k} = ''; % not found
                end
            else
                Nomatch.enzymes{k} = '';
            end
        end

        [foundIds, MWs, seqs, foundIdx] = GetUniProtMWSeqBatch(Nomatch.enzymes, 'PauseBetweenRequests', 0.2);
        matchedGenes_fromNomatch = Nomatch.Genes(foundIdx);
        enzymeConstraints.genes     = [enzymeConstraints.genes; matchedGenes_fromNomatch];
        enzymeConstraints.enzymes   = [enzymeConstraints.enzymes; foundIds];
        enzymeConstraints.mw        = [enzymeConstraints.mw; MWs];
        enzymeConstraints.sequence  = [enzymeConstraints.sequence; seqs];
    end

    %%
    % Step 7: parse gene-related rxns in the enzymeConstraints structure
    if strcmp(ecModeltype, 'integrated') || strcmp(ecModeltype, 'isozyme')
        enzymeConstraints.rxnEnzMat = zeros(numel(rxnWithGene),numel(enzymeConstraints.genes)); 
        for r=1:numel(rxnWithGene)
            rxnGenes   = model.genes(find(model.rxnGeneMat(rxnWithGene(r),:)));
            [~,locEnz] = ismember(rxnGenes,enzymeConstraints.genes); % Could also parse directly from rxnGeneMat, but some genes might be missing from Uniprot DB
            if locEnz ~= 0
                enzymeConstraints.rxnEnzMat(r,locEnz) = 1; %Assume 1 copy per subunit or enzyme, can be modified later
            end
        end
    elseif strcmp(ecModeltype, 'basic')
        enzymeConstraints.rxnEnzMat = zeros(numel(rxnWithGene),numel(enzymeConstraints.genes));
        for i=1:numRxnsOriginal
            index = rxnWithGene(i);
            gpr =model.grRules{index};
            gpr =strrep(gpr,'(','');
            gpr =strrep(gpr,')','');
            gpr =strrep(gpr,' or ',';');
            if (orCounts(i) == 0)
                geneNames = {gpr};
            else
                %Split the string into gene names
                geneNames=regexp(geneString,';','split');
            end
            %Loop through the isozymes and set the rxnGeneMat
            for j = 1:length(geneNames)
                %Find the gene in the gene list If ' and ' relationship, first
                %split the genes
                fnd = strfind(geneNames{j},' and ');
                if ~isempty(fnd)
                    andGenes=regexp(geneNames{j},' and ','split');
                    enzymeConstraints.rxnEnzMat(nextIndex,ismember(enzymeConstraints.genes,andGenes)) = 1; %should be subunit stoichoimetry
                else
                    enzymeConstraints.rxnEnzMat(nextIndex,ismember(enzymeConstraints.genes,geneNames(j)))=1;%should be subunit stoichoimetry
                end
                nextIndex = nextIndex + 1;
            end
        end
    end

    % Step 8: Add proteins as pseudometabolites
    if strcmp(ecModeltype, 'integrated')
        [proteinMets.mets, uniprotSortId] = unique(enzymeConstraints.enzymes);
        proteinMets.mets         = strcat('prot_',proteinMets.mets);
        proteinMets.metNames     = proteinMets.mets;
        proteinMets.compartments = 'c';
        if isfield(model,'metMiriams')
            proteinMets.metMiriams   = repmat({struct('name',{{'sbo'}},'value',{{'SBO:0000252'}})},numel(proteinMets.mets),1);
        end
        if isfield(model,'metCharges')
            proteinMets.metCharges   = zeros(numel(proteinMets.mets),1);
        end
        proteinMets.metNotes     = repmat({'Enzyme-usage pseudometabolite'},numel(proteinMets.mets),1);
        model = addMets(model,proteinMets);
    end
    
    % Step 9: Add protein pool pseudometabolite
    pool.mets         = 'prot_pool';
    pool.metNames     = pool.mets;
    pool.compartments = 'c';
    pool.metNotes     = 'Enzyme-usage protein pool';
    model = addMets(model,pool);
    
    % Step 10: Add protein usage reactions.
    if strcmp(ecModeltype, 'integrated')
        usageRxns.rxns            = strcat('usage_',proteinMets.mets);
        usageRxns.rxnNames        = usageRxns.rxns;
        usageRxns.mets            = cell(numel(usageRxns.rxns),1);
        usageRxns.stoichCoeffs    = cell(numel(usageRxns.rxns),1);
        for i=1:numel(usageRxns.mets)
            usageRxns.mets{i}         = {proteinMets.mets{i}, 'prot_pool'};
            usageRxns.stoichCoeffs{i} = [-1,1];
        end
        usageRxns.lb              = zeros(numel(usageRxns.rxns),1) - 1000;
        usageRxns.ub              = zeros(numel(usageRxns.rxns),1);
        usageRxns.rev             = ones(numel(usageRxns.rxns),1);
        usageRxns.grRules         = enzymeConstraints.genes(uniprotSortId);
        model = addRxns(model,usageRxns);
    end
    
    % Step 11: Add protein pool reaction (with open UB)
    poolRxn.rxns            = 'prot_pool_exchange';
    poolRxn.rxnNames        = poolRxn.rxns;
    poolRxn.mets            = {'prot_pool'};
    poolRxn.stoichCoeffs    = {-1};
    poolRxn.lb              = -1000;
    poolRxn.ub              = 0;
    poolRxn.rev             = 1;
    model = addRxns(model,poolRxn);
    
    fields = fieldnames(enzymeConstraints);
    for i = 1:length(fields)
        if isfield(model.enzymeConstraints, fields{i})
            model.enzymeConstraints.(fields{i}) = enzymeConstraints.(fields{i});
        end
    end
    ecModel = model;
    ecModel.enzymeConstraints.ecModeltype = ecModeltype;
end
