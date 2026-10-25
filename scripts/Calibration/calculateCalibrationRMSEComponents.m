
function componentRMSE = calculateCalibrationRMSEComponents( ...
    ecModel, bioRxn, cSource, basePath, org_name, ...
    modeltype, rmseMode, parameters)
% calculateCalibrationRMSEComponents
%   Calculates separate RMSE values for growth and exchange fluxes
%   using constrained growth experiments.
%
% Inputs:
%   ecModel   - Enzyme-constrained metabolic model
%   bioRxn    - Biomass reaction identifier
%   cSource   - Carbon source exchange reaction identifier
%   basePath  - Directory containing experimental data
%   org_name  - Organism identifier
%   modeltype - Model type (e.g., Bayesian, GECKO, or PRESTO)
%   rmseMode  - RMSE evaluation mode: 'flux' or 'carbon'
%   parameters - Optional parameter structure
%
% Output:
%   componentRMSE.rmseGrowth       - Growth-rate RMSE
%   componentRMSE.rmseExchangeFlux - Exchange-flux RMSE, excluding
%                                    the experimentally fixed substrate
%                                    uptake rate.
%
% Experimental and simulated values returned by abc_max are raw values.
% Carbon weighting is applied here when rmseMode is 'carbon'.

    % Initialize optional parameters
    if nargin < 8 || isempty(parameters)
        parameters = struct();
    end

    if nargin < 7 || isempty(rmseMode)
        rmseMode = 'flux';
    end

    if nargin < 6 || isempty(modeltype)
        modeltype = 'Bayesian';
    end

    % Retrieve missing inputs from the parameter structure
    if nargin < 5 || isempty(org_name)
        org_name = getRequiredParameter(parameters, 'org_name');
    end

    if nargin < 4 || isempty(basePath)
        basePath = getRequiredParameter(parameters, 'dataDir');
    end

    if nargin < 3 || isempty(cSource)
        cSource = getRequiredParameter(parameters, 'c_source');
    end

    if nargin < 2 || isempty(bioRxn)
        bioRxn = getRequiredParameter(parameters, 'bioRxn');
    end

    rmseMode = validatestring(rmseMode, {'flux', 'carbon'});

    % Add carbon-number information
    ecModel = addCarbonNum(ecModel, bioRxn);

    % Load constrained growth data
    growthData = readtable( ...
        fullfile(basePath, 'BayesianGrowthRates.tsv'), ...
        'FileType', 'text', ...
        'ReadRowNames', true);

    % Obtain experimental and simulated values
    % Column 1: substrate uptake
    % Column 2: growth rate
    % Columns 3-9: other exchange fluxes
    if ismember(string(modeltype), ["ECMpy", "GECKO", "PRESTO"])

        [~, exp_out, sim_out] = evaluateKcatRMSE( ...
            ecModel, growthData, [], [], ...
            bioRxn, cSource, [], org_name, [], rmseMode);

    else

        kcats = ecModel.enzymeConstraints.kcat;

        [~, exp_out, sim_out] = abc_max( ...
            ecModel, kcats, growthData, [], [], ...
            1, 1, 1, bioRxn, cSource, [], ...
            org_name, [], rmseMode);

    end

    % Handle unsuccessful simulations
    if isempty(sim_out) || all(~isfinite(sim_out(:)))

        componentRMSE.rmseGrowth = 999;
        componentRMSE.rmseExchangeFlux = 999;
        return;

    end

    % Remove experiments with invalid growth predictions
    validRows = isfinite(sim_out(:, 2));

    exp_out = exp_out(validRows, :);
    sim_out = sim_out(validRows, :);
    growthData = growthData(validRows, :);

    if isempty(sim_out)

        componentRMSE.rmseGrowth = 999;
        componentRMSE.rmseExchangeFlux = 999;
        return;

    end

    % Apply carbon weighting when required
    if strcmp(rmseMode, 'carbon')

        carbonNumbers = ecModel.excarbon(:);
        n = size(exp_out, 1);

        % Map substrate exchange reactions
        substrateNames = string(growthData{:, 1});
        [~, substrateIdx] = ismember( ...
            cellstr(substrateNames), ecModel.rxns);

        % Identify the biomass reaction
        biomassIdx = find(strcmp(ecModel.rxns, bioRxn), 1);

        % Identify the seven exchange reactions
        exchangeNames = growthData.Properties.VariableNames(4:10);
        [~, exchangeIdx] = ismember(exchangeNames, ecModel.rxns);

        % Construct carbon weights in the same column order as exp_out
        substrateWeights = carbonNumbers(substrateIdx(:));
        growthWeight = carbonNumbers(biomassIdx);
        exchangeWeights = reshape(carbonNumbers(exchangeIdx), 1, []);

        carbonWeights = [ ...
            substrateWeights, ...
            repmat(growthWeight, n, 1), ...
            repmat(exchangeWeights, n, 1)];

        % Match the zero-weight handling in rmsecal
        carbonWeights(carbonWeights == 0) = 1;

        exp_eval = exp_out .* carbonWeights;
        sim_eval = sim_out .* carbonWeights;

    else

        % Use raw, unweighted values
        exp_eval = exp_out;
        sim_eval = sim_out;

    end

    % Calculate growth-rate RMSE
    componentRMSE.rmseGrowth = calculateRMSE( ...
        exp_eval(:, 2), sim_eval(:, 2));

    % Calculate exchange-flux RMSE
    % Exclude column 1 because substrate uptake is experimentally fixed.
    % Columns 3-9 contain the remaining exchange fluxes.
    componentRMSE.rmseExchangeFlux = calculateRMSE( ...
        exp_eval(:, 3:9), sim_eval(:, 3:9));

end


function value = calculateRMSE(expValues, simValues)
% calculateRMSE
%   Calculates RMSE using finite experimental-simulation pairs.

    valid = isfinite(expValues) & isfinite(simValues);

    if ~any(valid(:))
        value = NaN;
        return;
    end

    residuals = expValues(valid) - simValues(valid);
    value = sqrt(mean(residuals.^2));

end


function value = getRequiredParameter(parameters, fieldName)
% getRequiredParameter
%   Retrieves a required field from the parameter structure.

    if isfield(parameters, fieldName) && ...
            ~isempty(parameters.(fieldName))
        value = parameters.(fieldName);
    else
        error('Missing required input: %s.', fieldName);
    end

end
