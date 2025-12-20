function model = loadJSONmodel(model, filename)
    % Load JSON file and decode its content
    raw = fileread(filename);
    data = jsondecode(raw);

    % Metabolic Network Components
    % Extract reaction IDs, metabolite IDs, stoichiometric matrix, bounds, etc.
    model.rxns = cellfun(@(x) x.id, data.reactions, 'UniformOutput', false);
    model.mets = {data.metabolites.id}';
    model.S = buildStoichiometricMatrix(model.mets, data.reactions);
    model.lb = cellfun(@(x) x.lower_bound, data.reactions, 'UniformOutput', false);
    model.ub = cellfun(@(x) x.upper_bound, data.reactions, 'UniformOutput', false);
    model.rev = zeros(length(model.rxns), 1);
    model.c = zeros(numel(model.rxns),1);  % JSON file does not contain the objective function info, so reset it.
    model.b = zeros(length(model.mets), 1);
    
    % Compartment Information
    % Get compartment IDs and names
    model.comps = fieldnames(data.compartments);
    model.compNames = cellfun(@(x) data.compartments.(x), model.comps, 'UniformOutput', false);

    % Gene Information
    % Extract gene IDs and build gene annotation structures
    model.genes = cellfun(@(x) x.id, data.genes, 'UniformOutput', false);
    model.geneMiriams = buildGeneMiriams(data.genes);
    model.geneShortNames = cellfun(@(x) x.name, data.genes, 'UniformOutput', false);
    
    % Reaction Annotations
    % Get reaction names, gene-reaction rules, and build the reaction-gene matrix
    model.rxnNames = cellfun(@(x) x.name, data.reactions, 'UniformOutput', false);
    model.grRules = cellfun(@(x) x.gene_reaction_rule, data.reactions, 'UniformOutput', false);
    model.rxnGeneMat = buildRxnGeneMat(model);
    % model.subSystems = {};  % To be determined
    model.eccodes = [];
    % model.eccodes = buildEccodes(data.reactions);
    model.rxnMiriams = buildRxnMiriams(data.reactions);
    
    % Metabolite Annotations
    % Extract metabolite names, compartments, formulas, annotations, and charges
    model.metNames = {data.metabolites.name}';
    model.metComps = {data.metabolites.compartment}';
    model.metFormulas = {data.metabolites.formula}';
    metAnnotation = {data.metabolites.annotation}';
    model.metMiriams = buildMetMiriams(metAnnotation);
    model.metCharges = {data.metabolites.charge}';

    % Enzyme Constraints
    % Initialize arrays for enzyme-constrained reactions
    enzyme_rxns = {};   % Reaction IDs with valid kcat values
    kcat_values = [];   % Corresponding kcat values
    rxnIdx = [];        % Indices of enzyme_rxns in model.rxns
    
    % Loop over each reaction in the JSON data
    for i = 1:length(data.reactions)
        currentKcat = data.reactions{i}.kcat;  % Get current reaction's kcat value
        if isnumeric(currentKcat) && ~isempty(currentKcat)
            % If kcat is numeric and not empty, record reaction ID and kcat
            enzyme_rxns{end+1} = data.reactions{i}.id;
            kcat_values(end+1,1) = currentKcat;
            [found, idx] = ismember(data.reactions{i}.id, model.rxns);
            if found
                rxnIdx(end+1,1) = idx;  % Record the index in model.rxns
            else
                rxnIdx(end+1,1) = NaN;  % Record NaN if not found
            end
        end
    end
    
    % Filter out NaN indices to ensure valid indexing
    valid = ~isnan(rxnIdx);
    rxnIdx = rxnIdx(valid);
    enzyme_rxns = enzyme_rxns(valid);
    kcat_values = kcat_values(valid);
    
    % Initialize additional arrays for enzyme constraints
    eccodes_values = cell(length(enzyme_rxns), 1);  % To store corresponding EC codes
    sourse_values = repmat({'ECMpy'}, length(enzyme_rxns), 1);  % kcat source filled with 'ECMpy'
    note_values = repmat({''}, length(enzyme_rxns), 1);           % Note field filled with empty strings
    
    % Build reaction-enzyme association matrix using gene-reaction rules
    % Use model.genes (not model.enzymeConstraints.genes as it is not created yet)
    rxnEnzMat = buildRxnGeneMat(model);
    
    % Extract EC codes corresponding to each enzyme-constrained reaction
    %{
    for i = 1:length(enzyme_rxns)
        eccodes_values{i} = model.eccodes{rxnIdx(i)};
    end
    %}

    % Preallocate arrays for enzyme concentrations and molecular weights
    concs_values = nan(length(model.genes), 1);
    mw_values = nan(length(model.genes), 1);
    % Create a cell array for enzyme sequences, one per gene, initialized as empty strings
    sequence_values = repmat({''}, length(model.genes), 1);
    
    % Retrieve enzyme UniProt IDs from geneMiriams if available
    enzymes_values = repmat({''}, length(model.genes), 1);
    for j = 1:length(model.genes)
        geneMiriam = model.geneMiriams{j};  % Get MIRIAM annotation for the current gene
        if isempty(geneMiriam)
            enzymes_values{j} = '';
            continue;
        end
        idx = find(contains(lower(geneMiriam.name), 'uniprot'), 1, 'first');
        if ~isempty(idx)
             enzymes_values{j} = geneMiriam.value{idx};
        end
    end
    
    % 注释信息
    % fieldNames = fieldnames(data.enzyme_constraint.kcat_MW);
    % fieldValues = cell(size(fieldNames));
    %{
    for i = 1:numel(fieldNames)
        fieldValues{i} = data.enzyme_constraint.kcat_MW.(fieldNames{i});
    end
    %}
    
    % Construct the enzyme constraints structure (1x1 struct)
    model.enzymeConstraints = struct('ecModeltype', 'simple',...
                                     'Ptot', data.enzyme_constraint.total_protein_fraction, ...             % Total protein content [g/gDW] (double)
                                     'f', data.enzyme_constraint.enzyme_mass_fraction, ...                  % Mass fraction of enzymes (double)
                                     'sigma', data.enzyme_constraint.average_saturation, ...                % Enzyme saturation factor (double)
                                     'rxns', {enzyme_rxns'}, ...                                             % Reactions involving enzymes
                                     'kcat', kcat_values/3600, ...                                               % Enzyme kinetic parameters
                                     'source', {sourse_values}, ...                                         % kcat source
                                     'notes', {note_values}, ...            
                                     'eccodes', {eccodes_values}, ...                                       % Enzyme ec number
                                     'concs', concs_values, ...                                             % Enzyme concentration
                                     'genes', {model.genes}, ...                                            % Enzyme coding gene
                                     'enzymes', {enzymes_values}, ...                                       % Enzyme uniprot id
                                     'mw', mw_values, ...                                                   % Enzyme relative molecular weight
                                     'sequence', {sequence_values}, ...                                                    % Enzyme sequences
                                     'rxnEnzMat', rxnEnzMat);          
end