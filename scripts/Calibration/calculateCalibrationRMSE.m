function rmseResult = calculateCalibrationRMSE(ecModel, useConstraint, useUnconstrained, useC13Flux, modeltype, parameters)
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

    if nargin < 6 || isempty(parameters), parameters = ParameterManager.getParams(); end

    bioRxn = parameters.bioRxn; cSource = parameters.c_source;
    basePath = parameters.dataDir; org_name = parameters.org_name;
    kcats = ecModel.enzymeConstraints.kcat;

    % [Fix] Pre-calculate Carbon Numbers
    if ~isfield(ecModel, 'excarbon')
        ecModel = addCarbonNum(ecModel, bioRxn);
    end

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
        [RMSE_con, ~, ~]    = evaluateKcatRMSE(ecModel, growthData, [], [], bioRxn, cSource, [], org_name, C13ReactionMap);
        [RMSE_uncon, ~, ~]  = evaluateKcatRMSE(ecModel, [], unconsData, [], bioRxn, cSource, [], org_name, C13ReactionMap);
        [RMSE_13C, ~, ~]    = evaluateKcatRMSE(ecModel, [], [], C13Fluxdata, bioRxn, cSource, [], org_name, C13ReactionMap);
        [RMSE_final, ~, ~]  = evaluateKcatRMSE(ecModel, growthData, unconsData, C13Fluxdata, bioRxn, cSource, [], org_name, C13ReactionMap);
    else
        [RMSE_con, ~, ~]    = abc_max(ecModel, kcats, growthData, [], [], 1, 1, 1, bioRxn, cSource, [], org_name, C13ReactionMap);
        [RMSE_uncon, ~, ~]  = abc_max(ecModel, kcats, [], unconsData, [], 1, 1, 1, bioRxn, cSource, [], org_name, C13ReactionMap);
        [RMSE_13C, ~, ~]    = abc_max(ecModel, kcats, [], [], C13Fluxdata, 1, 1, 1, bioRxn, cSource, [], org_name, C13ReactionMap);
        [RMSE_final, ~, ~]  = abc_max(ecModel, kcats, growthData, unconsData, C13Fluxdata, 1, 1, 1, bioRxn, cSource, [], org_name, C13ReactionMap);
    end

    rmseResult.rmseConstraint = RMSE_con;
    rmseResult.rmseUnconstrained = RMSE_uncon;
    rmseResult.rmse13C = RMSE_13C;
    rmseResult.rmseTotal = RMSE_final;
end





