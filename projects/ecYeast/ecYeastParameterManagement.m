function obj = ecYeastParameterManagement()
% TEMPLATE User parameter template file (Template.m)
%   Description: Modify your project parameters below, then retrieve the structure using:
%       obj = Template();
%   The configuration will be available in obj.params

    % ===== Project paths =====
    obj.params.path             = fullfile('C:\chenqilin\projects\ECOMAP', 'projects', 'ecYeast');
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
    obj.params.projectName = 'ecYeast';        % Project name
    obj.params.InitialModel = 'Yeast9.yml';    % Initial model
    obj.params.modeltype  = 'Tradition';        % Model type: 'ECOMAP', 'Tradition', 'sMOMENT', 'ECMpy', 'GECKO'

    % ===== Enzyme constraint parameters =====
    % Adjust based on the target organism and literature
    obj.params.sigma   = 0.5;  % Average enzyme saturation
    obj.params.Ptot    = 0.5;  % Total protein content [g/gDCW]
    obj.params.f       = 0.5;  % Enzyme fraction of total protein [g/g protein]

    % ===== Organism metadata =====
    obj.params.org_name        = 'Saccharomyces cerevisiae';
    obj.params.uniprot.type        = 'proteome';                % 'taxonomy' or 'proteome'
    obj.params.uniprot.ID          = 'UP000002311';             % Uniprot taxonomy ID
    obj.params.uniprot.geneIDfield = 'gene_oln';                % Gene ID field
    obj.params.uniprot.reviewed    = true;                      % Use only reviewed entries if true
    obj.params.taxonomicID         = '559292';                  % Get complex data

    % ===== Core reaction IDs =====
    obj.params.c_source = 'r_1714';                           % Carbon source exchange reaction
    obj.params.bioRxn   = 'r_2111';                           % Biomass reaction ID

    % ===== PRESTO OPTION  =====
    % Specify whether the approach should be run parallelized

    obj.params.PRESTO.runParallel = true;
    obj.params.PRESTO.ncpu = 8;
    obj.params.PRESTO.nIter = 50;                               % set the number of iterations of k-fold cross-validation
    % correction parameters
    obj.params.PRESTO.epsilon = 1e5;
    obj.params.PRESTO.lambda = 1e-5;
    obj.params.PRESTO.theta = 0.6;
    obj.params.PRESTO.GAM = repelem(75.5522,31);
    obj.params.PRESTO.comptype = {'P'};
    obj.params.PRESTO.comps    = {'ala__L_c'   89.09      'P'     % A     Alanine         ala
                                  'cys__L_c'  121.16      'P'     % C     Cysteine        cys
                                  'asp__L_c'  133.11      'P'     % D     Aspartic acid   asp
                                  'glu__L_c'  147.13      'P'     % E     Glutamic acid   glu
                                  'phe__L_c'  165.19      'P'     % F     Phenylalanine   phe
                                  'gly_c'      75.07      'P'     % G     Glycine         gly
                                  'his__L_c'  155.15      'P'     % H     Histidine       his
                                  'ile__L_c'  131.17      'P'     % I     Isoleucine      ile
                                  'lys__L_c'  146.19      'P'     % K     Lysine          lys
                                  'leu__L_c'  131.17      'P'     % L     Leucine         leu
                                  'met__L_c'  149.21      'P'     % M     Methionine      met
                                  'asn__L_c'  132.12      'P'     % N     Asparagine      asn
                                  'pro__L_c'  115.13      'P'     % P     Proline         pro
                                  'gln__L_c'  146.14      'P'     % Q     Glutamine       gln
                                  'arg__L_c'  174.2       'P'     % R     Arginine        arg
                                  'ser__L_c'  105.09      'P'     % S     Serine          ser
                                  'thr__L_c'  119.12      'P'     % T     Threonine       thr
                                  'val__L_c'  117.15      'P'     % V     Valine          val
                                  'trp__L_c'  204.23      'P'     % W     Tryptophan      trp
                                  'tyr__L_c'  181.19      'P'     % Y     Tyrosine        tyr
                                                              };  
    obj.params.PRESTO.pol_cost = [37.7,12.8,26,26];
end
