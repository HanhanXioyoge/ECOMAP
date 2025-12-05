%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%  
%%% Configuration file for PRESTO                                       %%%
%%%  -execute this script in the top-level directory of the approach    %%%
%%%                                                                     %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify organism name
orgName = 'Escherichia coli';
%Model basename used
orgBasename = 'eciML1515';
% Growth-associated maintenance (GAM); set to NaN if it should be fitted
GAM = repelem(75.5522,31);
% mass fraction of all proteins included in the model (see GECKO documentation)
f = 0.52102;
% mass fraction of unmeasured proteins according to PAX DB (also discounting
% proteins that do not have a measured abundance across all conditions)
f_n = NaN;
% average in vitro enzyme saturation (fitted in GECKO)
sigma = 0.5;
% correction factor for protein abudances
protCorrFact = NaN;

% solver for linear optimization
cobraSolver = 'gurobi';

% prefixes for enzyme metabolites and enzyme usage reaction in the GECKO
% model
enzMetPfx = 'prot_';
enzRxnPfx = 'prot_';
% Specify whether the approach should be run parallelized
runParallel = true;
ncpu = 8;
% set the number of iterations of k-fold cross-validation
nIter = 50;

mwFile = fullfile(findECOMAProot, 'scripts', 'database', 'max_MW.txt');

% correction parameters
epsilon = 100000;
lambda = 1e-05;
theta = 0.6; % relative error

% check if COBRA toolbox is installed
try
    changeCobraSolver(cobraSolver,'LP',0);
    changeCobraSolverParams('LP', 'feasTol', 1e-9);
catch
    error('COBRA toolbox is not installed')
end

% get maximum kcat value
maxKcatFile = fullfile(findECOMAProot, 'scripts', 'database', 'max_KCAT.txt');
K = retrieveMaxKcat(maxKcatFile,orgName);
clear maxKcatFile

