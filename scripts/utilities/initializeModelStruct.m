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
                                     'PDB', [], ...             % Enzyme structure
                                     'rxnEnzMat', []);          
end