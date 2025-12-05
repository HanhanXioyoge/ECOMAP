function startPRESTO(model, parameters)

    % -------- Parameters --------
    if nargin < 2 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end

    if ~isfield(model,'enzymeConstraints') || isempty(model.enzymeConstraints)
        error('The model does not contain the enzymeConstraints structure.');
    end
    
    dataDir  = parameters.dataDir;
    modelDir = parameters.modelDir;
    batchModelDir = fullfile(modelDir, 'PRESTO_batchModel');

    % -------- Create Directories if Not Exist --------
    if ~isfolder(batchModelDir)
        mkdir(batchModelDir);
    end

    % -------- Set Log Directory --------
    diary(fullfile(dataDir, [regexprep(model.name, '-', '') '_startPRESTO.log']));

    % -------- Prepare Conditions --------
    [condNames, E, expVal, nutrExch, P] = getconditions(dataDir, [], true);
    
    % Randomize conditions to avoid order bias
    rng(2021);  % Set random seed for reproducibility
    cond_idx = randsample(1:length(condNames), length(condNames));  % Randomly sample condition indices

    % Loop through each condition and enhance model
    for i = 1:length(cond_idx)
        tmp_idx = cond_idx(i);  % Get condition index

        % Update parameters for current condition
        load(fullfile(dataDir, 'parameters.mat'));
        parameters.Ptot = P(tmp_idx);
        parameters.gR_exp = expVal(tmp_idx);
        save(fullfile(dataDir, 'parameters.mat'), 'parameters');
        
        % Update chemostat data
        [match, idx] = ismember({'r_1714', 'r_1992', 'r_1672'}, nutrExch.Properties.RowNames);
        header = readtable(fullfile(dataDir, 'databases/chemostatData_temp.tsv'), 'FileType', 'text', 'Delimiter', '\t').Properties.VariableNames;
        chemostat_m = [expVal(tmp_idx), abs(table2array(nutrExch(idx, tmp_idx)))'];
        writetable(array2table(chemostat_m, 'VariableNames', header), fullfile(dataDir, 'databases/chemostatData.tsv'), 'FileType', 'text', 'Delimiter', '\t');
        
        % Enhance the GEM model with updated parameters
        enhanceGEM(model, 'COBRA', true, ['ecYeast_', condNames{tmp_idx}], '8');
    end

    % Restore original chemostat data
    if isfile(fullfile(dataDir, 'databases/chemostatData_temp.tsv'))
        movefile(fullfile(dataDir, 'databases/chemostatData_temp.tsv'), fullfile(dataDir, 'databases/chemostatData.tsv'));
    end

    % Close the log file
    diary off;
end












if ~isfolder('Logs')
    mkdir Logs
end
diary('Logs/ecoli_getcondec.log')
%Load GEM
load('Data/ecoli-GEM.mat')
%read in chen data
[condNames, E, expVal, nutrExch, P] = readDavidi2016([], topDir);
cd([geckoDir '/geckomat'])
load('../databases/parameters.mat')
 %check for correct organism 
 if ~strcmp(parameters.org_name, 'escherichia coli')
     error('GECKO parameters do not belong to Ecoli model, check if scripts have been copied from ecModels folder')
 end
  %save initial parameters
 sav_Ptot=parameters.Ptot;
 sav_gR_exp=parameters.gR_exp;
 clear parameters
%select conditions to build model from 
rng(2021)
cond_idx=randsample(1:length(condNames), length(condNames));
%For each selected condition
for i=start_i:length(cond_idx)
    tmp_idx=cond_idx(i)
    %Update the parameters object
    load('../databases/parameters.mat')
    parameters.Ptot=P(tmp_idx);
    parameters.gR_exp=expVal(tmp_idx);
    save('../databases/parameters.mat', 'parameters')
    clear parameters
    
    enhanceGEM(model, 'COBRA', unmod, ['ecEcoli_', condNames{tmp_idx}], '1')
end

 load('../databases/parameters.mat')
 parameters.Ptot=sav_Ptot;
 parameters.gR_exp=sav_gR_exp;
 clear parameters
cd(topDir)
diary off
end