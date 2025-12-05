%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%  
%%% Configuration file for PRESTO                                       %%%
%%%  -execute this script in the top-level directory of the approach    %%%
%%%                                                                     %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify organism name
orgName = '${org_name}';
%Model basename used
orgBasename = '${id}';
% Growth-associated maintenance (GAM); set to NaN if it should be fitted
GAM = NaN;
% mass fraction of all proteins included in the model (see GECKO documentation)
f = ${f};
% mass fraction of unmeasured proteins according to PAX DB (also discounting
% proteins that do not have a measured abundance across all conditions)
f_n = ;
% average in vitro enzyme saturation (fitted in GECKO)
sigma = ${sigma};
% correction factor for protein abudances
protCorrFact = NaN;

% solver for linear optimization
cobraSolver = 'gurobi';

% prefixes for enzyme metabolites and enzyme usage reaction in the GECKO
% model
enzMetPfx = 'prot_';
enzRxnPfx = 'prot_';
% Specify whether the approach should be run parallelized
runParallel = ${runParallel};
ncpu = ${ncpu};
% set the number of iterations of k-fold cross-validation
nIter = ${nIter};

mwFile = fullfile(findECOMAProot, 'scripts', 'database', 'max_MW.txt');

% correction parameters
epsilon = ${epsilon};
lambda = ${lambda};
theta = ${theta}; % relative error

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