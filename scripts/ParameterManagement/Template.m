function obj = KEY_Template()
% TEMPLATE User parameter template file (Template.m)
%   Description: Modify your project parameters below, then retrieve the structure using:
%       obj = Template();
%   The configuration will be available in obj.params

    % ===== Project root directory path =====
    obj.params.path      = fullfile('KEY_PATH', 'KEY_NAME');

    % ===== All model file =====
    % Supports JSON, SBML (.xml/.sbml), or other formats
    obj.params.modelDir  = fullfile(obj.params.path, 'models');

    % ===== Data directory =====
    obj.params.dataDir  = fullfile(obj.params.path, 'data');

    % ===== Output directory =====
    obj.params.outputDir = fullfile(obj.params.path, 'analysis');

    % ===== Project information =====
    obj.params.projectName = 'KEY_NAME';        % Project name
    obj.params.InitialModel = '';    % Initial model
    obj.params.modeltype  = '';        % Model type: 'ECOMAP', 'Tradition', 'sMOMENT', 'ECMpy', 'GECKO'

    % ===== Enzyme constraint parameters =====
    % Adjust based on the target organism and literature
    obj.params.sigma   = ;  % Average enzyme saturation
    obj.params.Ptot    = ;  % Total protein content [g/gDCW]
    obj.params.f       = ;  % Enzyme fraction of total protein [g/g protein]

    % ===== Organism metadata =====
    obj.params.org_name        = '';
    obj.params.uniprot.type        = '';                % 'taxonomy' or 'proteome'
    obj.params.uniprot.ID          = '';                % Uniprot taxonomy ID
    obj.params.uniprot.geneIDfield = '';                % Gene ID field
    obj.params.uniprot.reviewed    = true;              % Use only reviewed entries if true
    obj.params.taxonomicID         = '';                % Get complex data

    % ===== Core reaction IDs =====
    obj.params.c_source = '';                           % Carbon source exchange reaction
    obj.params.bioRxn   = '';                           % Biomass reaction ID

    % ===== PRESTO OPTION  =====
    % Specify whether the approach should be run parallelized
    obj.params.PRESTO.runParallel = ;                   % true/false
    obj.params.PRESTO.ncpu = ;
    obj.params.PRESTO.nIter = ;                                      % set the number of iterations of k-fold cross-validation
    % correction parameters
    obj.params.PRESTO.epsilon = ;
    obj.params.PRESTO.lambda = ;
    obj.params.PRESTO.theta = ;
end
