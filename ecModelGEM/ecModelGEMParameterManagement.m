function obj = ecModelGEMParameterManagement()
% TEMPLATE User parameter template file (Template.m)
%   Description: Modify your project parameters below, then retrieve the structure using:
%       obj = Template();
%   The configuration will be available in obj.params

    % ===== Project root directory path =====
    obj.params.path      = fullfile('D:\project\ECOMAP\ECOMAP', 'ecModelGEM');

    % ===== All model file =====
    % Supports JSON, SBML (.xml/.sbml), or other formats
    obj.params.modelDir  = fullfile(obj.params.path, 'models');

    % ===== Data directory =====
    obj.params.dataDir  = fullfile(obj.params.path, 'data');

    % ===== Output directory =====
    obj.params.outputDir = fullfile(obj.params.path, 'analysis');

    % ===== Project information =====
    obj.params.projectName = 'ecModelGEM';        % Project name
    obj.params.InitialModel = '';    % Initial model
    obj.params.modeltype  = 'Tradition';        % Model type: 'ECOMAP', 'Tradition', 'sMOMENT', 'ECMpy', 'GECKO'

    % ===== Enzyme constraint parameters =====
    % Adjust based on the target organism and literature
    obj.params.sigma   = 0.5;  % Average enzyme saturation
    obj.params.Ptot    = 0.55;  % Total protein content [g/gDCW]
    obj.params.f       = 0.55;  % Enzyme fraction of total protein [g/g protein]

    % ===== Organism metadata =====
    obj.params.org_name        = 'Escherichia coli';
    obj.params.uniprot.type        = 'proteome';                % 'taxonomy' or 'proteome'
    obj.params.uniprot.ID          = 'UP000000625';                % Uniprot taxonomy ID
    obj.params.uniprot.geneIDfield = 'gene_oln';                % Gene ID field
    obj.params.uniprot.reviewed    = true;              % Use only reviewed entries if true
    obj.params.taxonomicID         = '83333';                % Get complex data

    % ===== Core reaction IDs =====
    obj.params.c_source = 'EX_glc__D_e';                           % Carbon source exchange reaction
    obj.params.bioRxn   = 'biomass';                           % Biomass reaction ID

    % ===== PRESTO OPTION  =====
    % Specify whether the approach should be run parallelized
    obj.params.PRESTO.runParallel = false;                   % true/false
    obj.params.PRESTO.ncpu = 4;
    obj.params.PRESTO.nIter = 100;                                      % set the number of iterations of k-fold cross-validation
    % correction parameters
    obj.params.PRESTO.epsilon = 0.1;
    obj.params.PRESTO.lambda = 0.5;
    obj.params.PRESTO.theta = 1;
end
