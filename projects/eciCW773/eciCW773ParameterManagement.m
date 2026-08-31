function obj = eciCW773ParameterManagement()
% TEMPLATE User parameter template file (Template.m)
%   Description: Modify your project parameters below, then retrieve the structure using:
%       obj = Template();
%   The configuration will be available in obj.params

    % ===== Project paths =====
    obj.params.path             = fullfile('C:\chenqilin\projects\ECOMAP', 'projects', 'eciCW773');
    obj.params.projectDir       = obj.params.path;
    obj.params.projectJson      = fullfile(obj.params.projectDir, 'project.json');
    obj.params.parameterManager = fullfile(obj.params.projectDir, [mfilename '.m']);

    % ===== Project module directories =====
    obj.params.modelsDir         = fullfile(obj.params.projectDir, 'models');
    obj.params.reconstructionDir = fullfile(obj.params.projectDir, 'reconstruction');
    obj.params.calibrationDir    = fullfile(obj.params.projectDir, 'calibration');
    obj.params.analysisDir       = fullfile(obj.params.projectDir, 'analysis');
    obj.params.designDir         = fullfile(obj.params.projectDir, 'design');

    % ===== Project information =====
    obj.params.projectName = 'eciCW773';        % Project name
    obj.params.InitialModel = 'iCW773.xml';    % Initial model
    obj.params.modeltype  = 'Tradition';        % Model type: 'ECOMAP', 'Tradition', 'sMOMENT', 'ECMpy', 'GECKO'

    % ===== Enzyme constraint parameters =====
    % Adjust based on the target organism and literature
    obj.params.sigma   = 0.5;  % Average enzyme saturation
    obj.params.Ptot    = 0.56;  % Total protein content [g/gDCW]
    obj.params.f       = 0.46;  % Enzyme fraction of total protein [g/g protein]

    % ===== Organism metadata =====
    obj.params.org_name        = 'Corynebacterium glutamicum';
    obj.params.uniprot.type        = 'proteome';                % 'taxonomy' or 'proteome'
    obj.params.uniprot.ID          = 'UP000000582';                % Uniprot taxonomy ID
    obj.params.uniprot.geneIDfield = 'gene_oln';                    % Gene ID field
    obj.params.uniprot.reviewed    = false;                    % Use only reviewed entries if true
    obj.params.taxonomicID         = '196627';                % Get complex data

    % ===== Core reaction IDs =====
    obj.params.c_source = 'EX_glc_e';                           % Carbon source exchange reaction
    obj.params.bioRxn   = 'CG_biomass_cgl_ATCC13032';                           % Biomass reaction ID

    % ===== PRESTO OPTION  =====
    % Specify whether the approach should be run parallelized
    obj.params.PRESTO.runParallel = false;              % true/false
    obj.params.PRESTO.ncpu = 1;
    obj.params.PRESTO.nIter = 50;                       % set the number of iterations of k-fold cross-validation
    % correction parameters
    obj.params.PRESTO.epsilon = 1e5;
    obj.params.PRESTO.lambda = 1e-5;
    obj.params.PRESTO.theta = 0.6;
end
