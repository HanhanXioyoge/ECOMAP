function obj = ecHumanParameterManagement()
% TEMPLATE User parameter template file (Template.m)
%   Description: Modify your project parameters below, then retrieve the structure using:
%       obj = Template();
%   The configuration will be available in obj.params

    % ===== Project root directory path =====
    obj.params.path      = fullfile('D:\project\ECOMAP\ECOMAP', 'ecHuman');

    % ===== All model file =====
    % Supports JSON, SBML (.xml/.sbml), or other formats
    obj.params.modelDir  = fullfile(obj.params.path, 'models');

    % ===== Data directory =====
    obj.params.dataDir  = fullfile(obj.params.path, 'data');

    % ===== Output directory =====
    obj.params.outputDir = fullfile(obj.params.path, 'analysis');

    % ===== Project information =====
    obj.params.projectName = 'ecHuman';         % Project name
    obj.params.InitialModel = 'Human-GEM.xml';               % Initial model
    obj.params.modeltype  = 'Tradition';        % Model type: 'ECOMAP', 'Tradition', 'sMOMENT', 'ECMpy', 'GECKO'

    % ===== Enzyme constraint parameters =====
    % Adjust based on the target organism and literature
    obj.params.sigma   = 0.1;                   % Average enzyme saturation
    obj.params.Ptot    = 0.505717472;           % Total protein content [g/gDCW]
    obj.params.f       = 0.412;                 % Enzyme fraction of total protein [g/g protein]

    % ===== Organism metadata =====
    obj.params.org_name        = 'Homo sapiens';
    obj.params.uniprot.type        = 'proteome';                % 'taxonomy' or 'proteome'
    obj.params.uniprot.ID          = 'UP000005640';             % Uniprot taxonomy ID
    obj.params.uniprot.geneIDfield = 'gene_primary';            % Gene ID field
    obj.params.uniprot.reviewed    = true;              % Use only reviewed entries if true
    obj.params.taxonomicID         = '9606';                % Get complex data

    % ===== Core reaction IDs =====
    obj.params.c_source = 'MAR09034';                           % Carbon source exchange reaction
    obj.params.bioRxn   = 'MAR13082';                           % Biomass reaction ID

    % ===== Enzyme pseudo-metabolite compartment =====
    obj.params.enzyme_comp = 'cytoplasm';               % e.g., 'cytoplasm', 'periplasm'
end
