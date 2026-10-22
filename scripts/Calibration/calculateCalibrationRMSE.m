function rmseResult = calculateCalibrationRMSE(ecModel, useConstraint, useUnconstrained, useC13Flux, bioRxn, cSource, basePath, org_name, modeltype, rmseMode, parameters)
% calculateCalibrationRMSE
%   Calculates comprehensive RMSE between experimental and simulated calibration data.
%   Supports three data types: constraint growth data, unconstrained growth data, and 13C flux.
%
% Input:
%   ecModel   - Enzyme-constrained model (with kcat already applied via UpdateSmatrix)
%   varargin  - Optional name-value pairs:
%       'projectName' - Project name for data files (default: ecModel.id)
%       'projectPath' - Project root path (default: findECOMAProot)
%       'bioRxn' - Biomass reaction name (default: 'biomass')
%       'cSource' - Carbon source exchange reaction (default: 'EX_glc__D_e')
%       'org_name' - Organism name (default: 'ecoli')
%       'useConstraint' - Include constraint data (default: true)
%       'useUnconstrained' - Include unconstrained data (default: true)
%       'use13Cflux' - Include 13C flux data (default: true)
%       'verbose' - Print detailed messages (default: false)
%
% Output:
%   rmseResult - Structure containing:
%       .rmseConstraint - RMSE for constraint data (NaN if not used)
%       .rmseUnconstrained - RMSE for unconstrained growth data (NaN if not used)
%       .rmse13C - RMSE for 13C flux data (NaN if not used)
%       .rmseTotal - Total weighted RMSE
%       .nConstraint - Number of constraint experiments
%       .nUnconstrained - Number of unconstrained experiments
%       .n13C - Number of 13C flux data points
%       .details - Detailed results structure
%
% Example:

    if nargin < 11 || isempty(parameters)
        if nargin >= 10 && ...
                ~isempty(bioRxn) && ...
                ~isempty(cSource) && ...
                ~isempty(basePath) && ...
                ~isempty(org_name) && ...
                ~isempty(modeltype) && ...
                ~isempty(rmseMode)
            parameters = struct();
        else
            error('Please provide all required inputs or initialize ParameterManager.');
        end
    end

    if nargin < 10 || isempty(rmseMode)
        rmseMode = 'flux';
    end
    
    if nargin < 5 || isempty(bioRxn), bioRxn = parameters.bioRxn; end
    if nargin < 6 || isempty(cSource), cSource = parameters.c_source; end
    if nargin < 7 || isempty(basePath), basePath = parameters.dataDir; end
    if nargin < 8 || isempty(org_name), org_name = parameters.org_name; end
    if nargin < 9 || isempty(modeltype), modeltype = 'Bayesian'; end
    
    kcats = ecModel.enzymeConstraints.kcat;

    % Recalculate after model transformations so stale exchange weights are
    % not reused during evaluation.
    ecModel = addCarbonNum(ecModel, bioRxn);

    % Data Loading
    growthData = [];  if useConstraint, growthData = readtable(fullfile(basePath,'BayesianGrowthRates.tsv'), 'FileType', 'text', 'ReadRowNames', true); end
    unconsData = [];  if useUnconstrained, unconsData = readtable(fullfile(basePath,'UnconstrainedMaxGrowth.tsv'), 'FileType', 'text', 'ReadRowNames', true); end

    % 13C Flux Data - Use dedicated loading function
    C13Fluxdata = [];
    C13ReactionMap = [];
    if useC13Flux
        c13File = fullfile(basePath, '13CFluxdata.tsv');
        if exist(c13File, 'file')
            C13Fluxdata = load13CData(c13File, ecModel);
            % Build reaction mapping once for all conditions
            [C13ReactionMap, validationReport] = buildC13ReactionMap(C13Fluxdata, ecModel);
            fprintf('[bayesianTuning] 13C reaction mapping built: %d/%d reactions matched\n', ...
                validationReport.matchedCount, validationReport.totalReactions);
        else
            warning('13C flux data file not found: %s', c13File);
        end
    end
    if ismember(modeltype, ["ECMpy", "GECKO", "PRESTO"])
        % Use simplified RMSE evaluation for calibrated models
        [RMSE_con, ~, ~]    = evaluateKcatRMSE(ecModel, growthData, [], [], bioRxn, cSource, [], org_name, C13ReactionMap, rmseMode);
        [RMSE_uncon, ~, ~]  = evaluateKcatRMSE(ecModel, [], unconsData, [], bioRxn, cSource, [], org_name, C13ReactionMap, rmseMode);
        [RMSE_13C, ~, ~]    = evaluateKcatRMSE(ecModel, [], [], C13Fluxdata, bioRxn, cSource, [], org_name, C13ReactionMap, rmseMode);
        [RMSE_final, ~, ~]  = evaluateKcatRMSE(ecModel, growthData, unconsData, C13Fluxdata, bioRxn, cSource, [], org_name, C13ReactionMap, rmseMode);
    else
        [RMSE_con, ~, ~]    = abc_max(ecModel, kcats, growthData, [], [], 1, 1, 1, bioRxn, cSource, [], org_name, C13ReactionMap, rmseMode);
        [RMSE_uncon, ~, ~]  = abc_max(ecModel, kcats, [], unconsData, [], 1, 1, 1, bioRxn, cSource, [], org_name, C13ReactionMap, rmseMode);
        [RMSE_13C, ~, ~]    = abc_max(ecModel, kcats, [], [], C13Fluxdata, 1, 1, 1, bioRxn, cSource, [], org_name, C13ReactionMap, rmseMode);
        [RMSE_final, ~, ~]  = abc_max(ecModel, kcats, growthData, unconsData, C13Fluxdata, 1, 1, 1, bioRxn, cSource, [], org_name, C13ReactionMap, rmseMode);
    end

    rmseResult.rmseConstraint = RMSE_con;
    rmseResult.rmseUnconstrained = RMSE_uncon;
    rmseResult.rmse13C = RMSE_13C;
    rmseResult.rmseTotal = RMSE_final;
end
