function [ecModel, final_kcats, inital_kcats, rmse_history, kcat_history] = bayesianTuning(ecModel, useconsraint, useunconsraint, saveCurrentstate,maxIterations, proc, numPerGeneration, rejectnum, Initialvar, parameters)
    if nargin < 10 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end
    if nargin < 9
        Initialvar = 1;
    end
    if nargin < 8
        rejectnum = 0.2;
    end
    if nargin < 7
        numPerGeneration = 128;
    end
    if nargin < 6
        proc = parameters.PRESTO.ncpu;
    end
    if nargin < 5
        maxIterations = 150;
    end
    if nargin < 4
        saveCurrentstate = true;
    end
    if nargin < 3
        useconsraint = true;
    end
    if nargin < 2
        useunconsraint = true;
    end
    bioRxn = parameters.bioRxn;
    c_source = parameters.c_source;
    basePath = parameters.dataDir;
    org_name = parameters.org_name;
    if useconsraint
        growthdata = readtable(fullfile(basePath,'BayesianGrowthRates.tsv'), 'FileType', 'text', 'ReadRowNames', true);
    else
        growthdata = [];
    end
    if useunconsraint
        UnconstrainedMaxGrowth = readtable(fullfile(basePath,'UnconstrainedMaxGrowth.tsv'), 'FileType', 'text', 'ReadRowNames', true);
    else
        UnconstrainedMaxGrowth = [];
    end

    if isempty(growthdata) && isempty(UnconstrainedMaxGrowth)
        error('No data was read for calibration')
    end

    rxn2block = [];
    kcats = ecModel.enzymeConstraints.kcat;
    inital_kcats = kcats;
    kcat_var = ones(length(ecModel.enzymeConstraints.kcat),1)*Initialvar;
    
    % Initialize RMSE and kcat history
    rmse_history = [];  
    kcat_history = [];  
    sampledgeneration = 1;

    % Define the path where state will be saved
    analysisPath = parameters.outputDir;

    % Check if we have a saved state (for resuming after interruption)
    if saveCurrentstate
        if useconsraint && useunconsraint
            stateFile = fullfile(analysisPath, 'bayesianState_Total.mat');
        elseif ~useconsraint && useunconsraint
            stateFile = fullfile(analysisPath, 'bayesianState_UnconstrainedMaxGrowth.mat');
        elseif useconsraint && ~useunconsraint
            stateFile = fullfile(analysisPath, 'bayesianState_GrowthRates.mat');
        end
        if exist(stateFile, 'file')
            load(stateFile, 'rmse_history', 'kcat_history', 'theta_100', 'kcat_100', 'sampledgeneration');
            disp('Resuming from saved state...');
        else
            theta_100 = [];
            kcat_100 = [];
        end
    else
        theta_100 = [];
        kcat_100 = [];
    end

    % First test
    D = abc_max_test(ecModel, kcats, growthdata, UnconstrainedMaxGrowth, 1, 1, 1, bioRxn, c_source, rxn2block, org_name);
    D_100 = D;
    
    while D > rejectnum && D_100 > 0.25
        if sampledgeneration <= maxIterations
            disp(['Running ' num2str(sampledgeneration) ' of ' num2str(maxIterations) ': D = ' num2str(D) ' D_100 = ' num2str(D_100)])
            
            % Save the previous iteration's state before running the next
            old = theta_100;
            kcat_old_100 = kcat_100;

            if sampledgeneration == 1
                sample_generation = 160;
            else
                sample_generation = numPerGeneration;
            end
            
            % generate one maxIterations sample of kcats
            kcat_random_all = arrayfun(@getrSample, kcats, kcat_var, repmat(sample_generation, length(kcats), 1),'UniformOutput', false);
            kcat_random_all = cell2mat(kcat_random_all);
            
            % Simulate and measure RMSE
            parfor j = 1:proc
                rmse_final = abc_max_test(ecModel, kcat_random_all, growthdata, UnconstrainedMaxGrowth, proc, sample_generation, j, bioRxn, c_source, rxn2block, org_name);
                new_tmp{j} = rmse_final;
            end
            new = cell2mat(new_tmp);
            theta = [new, old];
            kcat = [kcat_random_all, kcat_old_100];
            
            % Initialize an empty set to store the best 100 after each step
            [~, D_idx] = sort(theta, 'ascend');
            theta_100 = theta(D_idx(1:100));
            D = abs(theta_100(100) - theta_100(1));
            D_100 = theta_100(100);
            kcat_100 = kcat(:, D_idx(1:100));
            
            % Record RMSE and kcat history
            rmse_history = [rmse_history, mean(new, 'omitnan')];
            kcat_history = [kcat_history, kcat_100(:, 1)];
            
            % Recalculate the sigma and mu
            ss = num2cell(kcat_100', 1);
            [a, b] = arrayfun(@updateprior, ss);
            kcats = a';
            kcat_var = b';
            
            % Save the current state after each iteration
            if saveCurrentstate
                save(stateFile, 'rmse_history', 'kcat_history', 'theta_100', 'kcat_100', 'sampledgeneration');
            end
            
            sampledgeneration = sampledgeneration + 1;
        else
            D = rejectnum;
            D_100 = D;
        end
    end

    final_kcats = kcats;
    ecModel.enzymeConstraints.kcat = kcats;
end

function r = getrSample(mu,sigma,step,method)
% getrSample
%   Samples random kcats from a distribution.
%
% Input:
%   mu              Mean of distribution (data is logged to get a normal distr)
%   sigma           Std deviation of the distribution
%   step            Number of kcats to sample
%   method          shape of distribution: 'normal' or 'uniform'. 
%                   (Optional, default is 'normal')
% Output:
%   r               The sampled kcats
%
if nargin < 4
    method = 'normal';
end
if mu == 0
    r = zeros(1,step);
elseif strcmp(method,'normal')
    mutmp = log10(mu);
    %sigmatmp = log10(sigma);
    sigmatmp = sigma;
    pd = makedist('normal','mu',mutmp,'sigma',sigmatmp);
    t = truncate(pd,-3,7);
    r = random(t,1,step);
    r = 10.^(r);
elseif strcmp(method,'uniform')
    mutmp = log10(mu);
    sigmatmp = sigma;
    pd = makedist('uniform','lower',mutmp-sigmatmp,'upper',mutmp + sigmatmp);
    t = truncate(pd,-3,7);
    r = random(t,1,step);
    r = 10.^(r);
end

r(r<0) = 0;

end

function [a,b] = updateprior(x)
% updateprior
%   Calculates a new distribution from the selected kcat values
%
% Input:
%   x               kcats
% Output:
%   a               Mean
%   b               Std dev

x = x{:};
x(x == 0) = []; %remove zeros - they cannot be handled in the log transform below
if length(x) == 0
    %we don't have much choice in the two first cases - just set sigma to 1 - same as in the initial prior
    a = 0;
    b = 1;
elseif length(x) == 1
    a = x;
    b = 1;
else
    MIN_SIGMA = 0.05;
    pd = fitdist(log10(x),'Normal');
    a = 10^(pd.mu);
    b = max(pd.sigma, MIN_SIGMA);
end

end