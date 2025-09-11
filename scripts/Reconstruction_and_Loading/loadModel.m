function model = loadModel(filename, modeltype)
% loadModel - Load metabolic models in various formats and convert to unified structure
% Inputs:
%   filename     Model file name (optional, default: file selection dialog)
%   modeltype    Model type ('ECOMAP','COBRA','sMOMENT','ECMpy','GECKO')
%   EnzymeData   Enzyme data dictionary containing kcat and molecular weight
% Output:
%   model        Unified metabolic model structure with standardized fields

% Handle input parameters
if nargin<1 || isempty(filename)
    [filename, pathName] = uigetfile({'*.xml;*.json;*.yml;*.mat','Model Files (*.xml,*.json,*.yml,*.mat)'},...
                                    'Select model file');
    if filename == 0
        error('Model file selection required');
    end
    filename = fullfile(pathName, filename);
end

if ~isfile(filename)
    error('File not found: %s', filename);
end

if nargin < 2 || isempty(modeltype)
    opts = {'ECOMAP', 'COBRA', 'sMOMENT', 'ECMpy', 'GECKO'};
    [selected, ok] = listdlg('ListString', opts, 'SelectionMode', 'single',...
                           'PromptString', 'Select model type:',...
                           'ListSize', [200, 150]);
    if ~ok, error('Model type selection required'); end
    modeltype = opts{selected};
end

modeltype = upper(modeltype);
validTypes = {'ECOMAP','COBRA','SMOMENT','ECMPY','GECKO'};
if ~ismember(modeltype, validTypes)
    error('Unsupported model type: %s', modeltype);
end

% Initialize unified structure
model = initializeModelStruct();

% Reading information based on model file type
[~, fileNameWithoutExtension, FileExtension] = fileparts(filename);

% Load model based on FileExtension
if strcmp(FileExtension, '.xml')
    model = loadXML(model, filename);
elseif strcmp(FileExtension, '.json')
    model = loadJSON(model, filename);
elseif strcmp(FileExtension, '.yml')
    model = loadYML(model, filename);
else
    error(['Cannot process files of type ' FileExtension]);
end

% Standardize model structure
model = standardizeModel(model, fileNameWithoutExtension, modeltype);
end

%% Subfunction: Initialize model structure 
function model = initializeModelStruct()
% Initialize model structure with annotated fields
model = [];
model.id = [];
model.name = [];

% Biological context information
model.information = struct('organism', [], ...          % Source organism name (string)
                           'taxonomicID', [], ...       % NCBI taxonomic identifier (string)
                           'uniprot_type', [], ...      % UniProt identifier system used (string)
                           'uniprot_id', []);            % UniProt protein accessions (cell array)

% Metabolic network components
model.rxns = {};       % Reaction IDs (cell array)
model.mets = {};       % Metabolite IDs (cell array)
model.S = [];          % Stoichiometric matrix (sparse matrix)
model.lb = [];         % Lower flux bounds (double vector)
model.ub = [];         % Upper flux bounds (double vector)
model.rev = [];        % Reversibility flags (0/1 vector)
model.c = [];          % Objective coefficients (double vector)
model.b = [];          % Metabolite balance constraints (double vector)

% Compartment information
model.comps = {};      % Compartment IDs (cell array)
model.compNames = {};  % Compartment names (cell array)

% Gene information
model.genes = {};         % Gene IDs (cell array)
model.geneMiriams = {};   % MIRIAM annotations for genes (cell array)
model.geneShortNames = {};% Gene abbreviations (cell array)

% Reaction annotations
model.rxnNames = {};      % Reaction names (cell array)
model.grRules = {};       % Gene-reaction rules (cell array)
model.rxnGeneMat = [];    % Gene-reaction associations (sparse matrix)
% model.subSystems = {};    % Metabolic subsystem classifications (cell array)
model.eccodes = {};       % EC numbers (cell array)
model.rxnMiriams = {};    % MIRIAM annotations for reactions (cell array)

% Metabolite annotations
model.metNames = {};      % Metabolite names (cell array)
model.metComps = [];      % Compartment indices for metabolites (double vector)
model.metFormulas = {};   % Chemical formulas (cell array)
model.metMiriams = {};    % MIRIAM annotations for metabolites (cell array)
model.metCharges = [];    % Charge states (double vector)

% Enzyme constraints
model.enzymeConstraints = struct('ecModeltype',[],...       % The type of ecModel includes two types: complex and simple
                                 'Ptot', [], ...            % Total protein content [g/gDW] (double)
                                 'f', [], ...               % Mass fraction of enzymes (double)
                                 'sigma', [], ...           % Enzyme saturation factor (double)
                                 'rxns', [], ...            % Reactions involving enzymes
                                 'kcat', [], ...            % Enzyme kinetic parameters
                                 'source', [], ...          % kcat source
                                 'notes', [], ...            
                                 'eccodes', [], ...         % Enzyme ec number
                                 'concs', [], ...           % Enzyme concentration
                                 'genes', [], ...           % Enzyme coding gene
                                 'enzymes', [], ...         % Enzyme uniprot id
                                 'mw', [], ...              % Enzyme relative molecular weight
                                 'sequence', [], ...        % Enzyme sequences
                                 'rxnEnzMat', []);          
end

%% Subfunctions for model loading
function model = loadJSON(model, filename)
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
    model.rxnGeneMat = buildRxnGeneMat(model.genes, model.grRules);
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
    rxnEnzMat = buildRxnGeneMat(model.genes, model.grRules(rxnIdx));
    rxnEnzMat = full(rxnEnzMat);  % Convert sparse matrix to full double matrix
    
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
                                     'kcat', kcat_values, ...                                               % Enzyme kinetic parameters
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

function model = loadXML(model, filename)
    % Load XML model
    try
        xml_model = importModel(filename);
    catch
        error('required for loading SBML models');
    end
    
    % Correspondence using names
    fields = fieldnames(xml_model);
    
    for i = 1:length(fields)
        field = fields{i};
        if isfield(model, field)
            model.(field) = xml_model.(field);
        end
    end
end

function model = loadYML(model, filename)
    % Load YML model from file
    try
        yml_model = readYAMLmodel(filename);
    catch
        error('Required function for loading YAML models is missing or failed.');
    end
    
    % Get all fieldnames from the yml_model
    fields = fieldnames(yml_model);
    
    % Loop over each field in the YAML model and copy it to model if that field already exists in model
    for i = 1:length(fields)
        field = fields{i};
        if isfield(model, field)
            model.(field) = yml_model.(field);
        end
    end
    
    % Check if the yml_model contains the 'ec' field for enzyme constraints
    if isfield(yml_model, 'ec')
        % Determining ecModel type
        if ~isfield(yml_model, 'ec') || ~isfield(yml_model.ec, 'geckoLight')
            error('The field yml_model.ec.geckoLight does not exist.');
        % Assign ecModeltype based on geckoLight value
        elseif yml_model.ec.geckoLight == true
            model.enzymeConstraints.ecModeltype = 'simple';
        elseif yml_model.ec.geckoLight == false
            model.enzymeConstraints.ecModeltype = 'complex';
        else
            error('yml_model.ec.geckoLight must be logical.');
        end

        % Define the allowed fields to be assigned from yml_model.ec
        allowed_fields = {'rxns', 'kcat', 'source', 'notes', 'eccodes', 'concs', 'genes', 'enzymes', 'mw', 'sequence', 'rxnEnzMat'};
        % Get the fieldnames from the yml_model.ec structure
        ec_fields = fieldnames(yml_model.ec);
        % Iterate over each field in yml_model.ec
        for i = 1:length(ec_fields)
            field = ec_fields{i};
            % If the current field is in the allowed list, assign it to model.enzymeConstraints
            if ismember(field, allowed_fields)
                model.enzymeConstraints.(field) = yml_model.ec.(field);
            end
        end
    end
end

%% Helper functions
function model = standardizeModel(model, fileNameWithoutExtension, modeltype)
    % Core model identification
    model.id = fileNameWithoutExtension;          % Unique model identifier (string)
    model.name = modeltype;                       % Human-readable model name (string)
end

function S = buildStoichiometricMatrix(mets, reactions)
% BUILDSTOICHIOMETRICMATRIX Construct a stoichiometric matrix 
% from the structure of the reactive metabolite
% input：
%   mets       - Metabolite ID list (array of cells)
%   reactions  - An array of reaction structures, each element 
%                contains a metabolite substructure
% output：
%   S          - sparse stoichiometric matrix（mets × reactions）

numMets = numel(mets);
numRxns = numel(reactions);

rowIdx = [];
colIdx = [];
values = [];

for rxnIdx = 1:numRxns
    metabolites = reactions{rxnIdx}.metabolites;
    
    metIDs = fieldnames(metabolites);
    stoich = struct2array(metabolites); 
    
    if ~isnumeric(stoich)
        error('The stoichiometric number of %d of the reaction should be of numerical type', rxnIdx);
    end
    
    [~, metPositions] = ismember(metIDs, mets);
    
    valid = metPositions ~= 0;
    metPositions = metPositions(valid);
    stoich = stoich(valid);
    
    rowIdx = [rowIdx; metPositions];
    colIdx = [colIdx; rxnIdx * ones(length(metPositions), 1)];
    values = [values; stoich(:)]; 
end

S = sparse(rowIdx, colIdx, values, numMets, numRxns);

% dimensionality verification
assert(size(S,1) == numMets, 'Inconsistent metabolite dimensions: expected %d actual %d', numMets, size(S,1));
assert(size(S,2) == numRxns, 'Inconsistent response dimensions: expected %d actual %d', numRxns, size(S,2));
end

function rxnGeneMat = buildRxnGeneMat(genes, grRules)
% BUILDRXNGENEMAT Constructs the reaction-gene association matrix from gene rules.
% Inputs:
%   genes   - cell array containing gene IDs.
%   grRules - cell array containing reaction gene rules.
% Output:
%   rxnGeneMat - logical sparse reaction-gene association matrix.

% Validate inputs
numRxns = numel(grRules);
numGenes = numel(genes);

% Initialize sparse matrix
rxnGeneMat = sparse(numRxns, numGenes);

% Regular expression pattern to match gene IDs (adaptable for different naming conventions)
genePattern = '(\<\w+\>)'; % Matches word tokens separated by word boundaries

for i = 1:numRxns
    ruleStr = grRules{i};
    
    % Skip empty rules
    if isempty(ruleStr)
        continue;
    end
    
    % Extract all unique gene IDs
    [geneTokens, ~] = regexp(ruleStr, genePattern, 'tokens', 'match');
    uniqueGenes = unique([geneTokens{:}]);
    
    % Filter out logical operators 'and' and 'or'
    uniqueGenes = uniqueGenes(~ismember(lower(uniqueGenes), {'and','or'}));
    
    % Find the indices of genes in the provided gene list
    [found, geneIdx] = ismember(uniqueGenes, genes);
    
    if any(~found)
        warning('Reaction %d contains unrecognized gene(s): %s', i, strjoin(uniqueGenes(~found), ', '));
        geneIdx = geneIdx(found);
    end
    
    % Update the association matrix
    if ~isempty(geneIdx)
        rxnGeneMat(i, geneIdx) = 1;
    end
end

% Convert to logical sparse matrix
% rxnGeneMat = logical(rxnGeneMat);

% Validate matrix dimensions
assert(size(rxnGeneMat,1) == numRxns, 'Number of reactions does not match.');
assert(size(rxnGeneMat,2) == numGenes, 'Number of genes does not match.');
end

function geneMiriams = buildGeneMiriams(genes)
% BUILDGENEMIRIAMS Constructs MIRIAM annotation structures for genes.
%
% Input:
%   genes - A cell array of gene structures. Each gene may have an 'annotation' field.
%
% Output:
%   geneMiriams - A cell array where each element is a scalar structure with two fields:
%                 'name' and 'value'. The 'name' field is a cell array of annotation field
%                 names, and the 'value' field is a cell array of the corresponding values.

% Determine the number of genes
numGenes = numel(genes);
geneMiriams = cell(numGenes, 1);

for i = 1:numGenes
    if isfield(genes{i}, 'annotation')
        % Extract annotation field names and their corresponding values
        annFields = fieldnames(genes{i}.annotation);
        annValues = struct2cell(genes{i}.annotation);
        % Build a scalar structure with fields 'name' and 'value'
        miriamStruct.name = annFields;
        miriamStruct.value = annValues;
    else
        % If there is no annotation field, assign empty cell arrays
        miriamStruct.name = {};
        miriamStruct.value = {};
    end
    geneMiriams{i} = miriamStruct;
end
end

function rxnMiriams = buildRxnMiriams(rxns)
% BUILDRXNMIRIAMS Constructs MIRIAM annotation structures for reactions.
%
% Input:
%   rxns - A cell array of reaction structures. Each reaction may have an
%          'annotation' field.
%
% Output:
%   rxnMiriams - A cell array where each element is a scalar structure with two fields:
%                'name' and 'value'. For annotation fields that contain multiple elements,
%                the field name appears repeatedly (one for each element) and the values are
%                flattened into the 'value' cell array.

numRxns = numel(rxns);
rxnMiriams = cell(numRxns, 1);

for i = 1:numRxns
    if isfield(rxns{i}, 'annotation')
        annFields = fieldnames(rxns{i}.annotation);   
        % Initialize new cell arrays for the expanded field names and values.
        newAnnFields = {};
        newAnnValues = {};
        % Loop over each annotation field.
        for j = 1:length(annFields)
            currentField = annFields{j};
            currentValue = rxns{i}.annotation.(currentField);
            if iscell(currentValue)
                % If the value is a cell array, replicate the field name for each element.
                newAnnFields = [newAnnFields; repmat({currentField}, numel(currentValue), 1)];
                % Flatten the cell array and append each element.
                newAnnValues = [newAnnValues; currentValue(:)];
            else
                % Otherwise, just add the field name and its value once.
                newAnnFields = [newAnnFields; {currentField}];
                newAnnValues = [newAnnValues; {currentValue}];
            end
        end
        newAnnFields = cellfun(@(s) strrep(s, '_', '.'), newAnnFields, 'UniformOutput', false);
        % Build a scalar structure with fields 'name' and 'value'.
        miriamStruct.name = newAnnFields;
        miriamStruct.value = newAnnValues;
    else
        % If there is no annotation field, assign empty cell arrays.
        miriamStruct.name = {};
        miriamStruct.value = {};
    end
    rxnMiriams{i} = miriamStruct;
end
end

function metMiriams = buildMetMiriams(metAnnotation)
% BUILDMETMIRIAMS Constructs MIRIAM annotation structures for metabolites.
%
% Input:
%   metAnnotation - A cell array of metabolite annotation structures.
%                   Each cell contains a structure with various annotation fields.
%
% Output:
%   metMiriams   - A cell array where each element is a scalar structure 
%                  with two fields:
%                  'name'  - a cell array of annotation field names (with underscores replaced by dots),
%                  'value' - a cell array of the corresponding annotation values.

numMets = numel(metAnnotation);
metMiriams = cell(numMets, 1);

for i = 1:numMets
    % Initialize an empty structure with fields 'name' and 'value'
    miriamStruct.name = {};
    miriamStruct.value = {};
    
    % Process if the current metabolite annotation is a structure.
    if isstruct(metAnnotation{i})
        % Extract all annotation field names.
        annFields = fieldnames(metAnnotation{i});
        % Initialize new cell arrays for expanded field names and values.
        newAnnFields = {};
        newAnnValues = {};
        % Loop over each annotation field.
        for j = 1:length(annFields)
            currentField = annFields{j};
            currentValue = metAnnotation{i}.(currentField);
            if iscell(currentValue)
                % If the value is a cell array, replicate the field name for each element.
                newAnnFields = [newAnnFields; repmat({currentField}, numel(currentValue), 1)];
                % Flatten the cell array and append each element.
                newAnnValues = [newAnnValues; currentValue(:)];
            else
                % Otherwise, add the field name and its value once.
                newAnnFields = [newAnnFields; {currentField}];
                newAnnValues = [newAnnValues; {currentValue}];
            end
        end
        % Replace underscores with dots in newAnnFields.
        newAnnFields = cellfun(@(s) strrep(s, '_', '.'), newAnnFields, 'UniformOutput', false);
        % Build the scalar structure with fields 'name' and 'value'.
        miriamStruct.name = newAnnFields;
        miriamStruct.value = newAnnValues;
    end
    metMiriams{i} = miriamStruct;
end
end