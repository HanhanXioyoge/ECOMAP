function [solution,models,relError,changeTab,LP] = PRESTOforECOMAP(varargin)

p = parseInput(varargin);

% Initialize the parameters
models = p.Results.models;
expVal = p.Results.expVal;
E = p.Results.E;
enzMetPfx = p.Results.enzMetPfx;
enzRxnPfx = p.Results.enzRxnPfx;
epsilon = p.Results.epsilon;
theta = p.Results.theta;
lambda = p.Results.lambda;
enzBlackList = p.Results.enzBlackList;
K = 3600*p.Results.K;
includeUM = p.Results.includeUM;
pCorrFactor = p.Results.pCorrFactor;
negCorrFlag = p.Results.negCorrFlag;

% lambda = 1e-5;
% number of conditions
NCOND = size(E,2);
if numel(expVal) ~= NCOND
    error('Number of provided growth rates does not match the number of conditions')
elseif numel(models) ~= NCOND
    error('Number of provided models does not match the number of conditions')
end

if isstruct(models)
    models = {models};
end

% if pCorrFactor is a scalar, transform to vector with as many entries as
% models
if isscalar(pCorrFactor)
    pCorrFactor = pCorrFactor .* ones(numel(models),1);
end

% find commmon indices using first model as representative
model = models{1};

% find enzyme indices in model metabolites
[enzMetIdx,enzRxnIdx,poolMetIdx,poolRxnIdx] = findEnzIdx(model,enzMetPfx,enzRxnPfx);

% get row and column dimensions of stoichiometric matrix
[RDIM,CDIM] = size(model.S);
% determine the number of deltas from the number of enzyme usage reactions
allRxnIdx = 1:CDIM;
NDELTA    = numel(enzRxnIdx);

% 生化反应 = 所有反应 - （有 Δ 的反应）-（你选择剔除的池反应）
bioRxnIdx    = setdiff(allRxnIdx, enzRxnIdx);
only_bioRxnIdx = setdiff(bioRxnIdx,poolRxnIdx);
nBiochemRxns = numel(bioRxnIdx);
only_nBiochemRxns = numel(only_bioRxnIdx);

% find biomass reaction and protein pool reaction indices in models
bioIdx = findBioIdx(models,only_bioRxnIdx, nBiochemRxns);
ProteinpoolRxnIdx = findproteinpoolIdx(models,poolRxnIdx, nBiochemRxns);   % column
ProteinpoolMetIdx = findproteinpoolIdx(models,poolMetIdx, RDIM);           % row

% get blacklist indices (enzymes, which should be excluded from correction)
blackListIdx = ismember(model.mets(enzMetIdx),enzBlackList);

idxExclude = ~all(E>0,2);

% store original and minimum kcat values per enzyme (identical across all models)
[~,minKcat] = getKcatFromEcModel(model, enzMetPfx, enzRxnPfx);
MW = model.enzymeConstraints.mw; % mg/mmol

% initialize matrix and rhs for the final linear optimization program
nz = NCOND * sum(sum(model.S~=0)) + sum(any(E,2)) - sum(blackListIdx);
lpMatrix = spalloc(NCOND*RDIM, NCOND*nBiochemRxns+NDELTA, nz);
% right hand side
lpRHS = zeros(size(lpMatrix,1),1);
% constraint sense (COBRA style)
tmpSense = model.csense;
tmpSense(enzMetIdx) = 'L';
lpSense = repmat(tmpSense,NCOND,1);
clear model tmpSense
% LP objective (minimize deltas)
deltaColStartIdx = NCOND*nBiochemRxns+1;
lpObj = zeros(size(lpMatrix,2),1);

% weigh deltas by lambda, divided by the number of kcats that can be corrected
lpObj(deltaColStartIdx:end) =  lambda / (sum(any(E,2)) - sum(blackListIdx));

for i=1:numel(models) 
    tmpModel = models{i};
    % tmpModel.S(poolMetIdx, poolRxnIdx)
    tmpModel.S(enzMetIdx,enzRxnIdx) = 0;
    for j = 1:numel(enzMetIdx)
        enzRow    = enzMetIdx(j);
        enzKcatIdx = tmpModel.S(enzRow,:) < 0;

        if any(enzKcatIdx)
            if ismember(enzRow, enzMetIdx(idxExclude))
                tmpModel.S(enzRow,:) = 0;
            else

                tmpModel.S(enzRow, enzKcatIdx) = MW(j);
                rowGlobal   = (i-1)*RDIM + enzRow;
                colDelta_e  = deltaColStartIdx + (j-1);
                lpMatrix(rowGlobal, colDelta_e) = -E(j,i);
            end
        end
    end
    clear enzKcatIdx

    tmpSMat = tmpModel.S(:, bioRxnIdx);
    colMask = ismember(bioRxnIdx, only_bioRxnIdx);
    tmpSMat(poolMetIdx, colMask) = sum(tmpSMat(enzMetIdx, colMask), 1);

    for j = 1:numel(only_bioRxnIdx)
        bioRxnGlob = only_bioRxnIdx(j);
        colTmp     = find(bioRxnIdx == bioRxnGlob, 1);
    
        if isempty(colTmp)
            continue;
        end
    
        enzKcatIdx = tmpSMat(enzMetIdx,colTmp) > 0;
        if any(enzKcatIdx)
            minkcat  = min(minKcat(enzKcatIdx));
            coeff    = -tmpSMat(poolMetIdx,colTmp) / minkcat;
            tmpSMat(poolMetIdx,colTmp) = coeff;
        end
    end

    rowIdx = (i-1)*RDIM + 1 : i*RDIM;
    colIdx = (i-1)*nBiochemRxns + 1 : i*nBiochemRxns;
    % Determined size(tmpSMat) = rowIdx*colIdx
    lpMatrix(rowIdx, colIdx) = tmpSMat;
    
    tmpRHS = tmpModel.b;
    tmpRHS(enzMetIdx) = minKcat .* E(:,i);
    lpRHS(rowIdx) = tmpRHS;
end
lpMatrix = sparse(lpMatrix);
clear rowIdx colIdx enzKcatIdx tmpSMat tmpDeltaMat tmpModel tmpRHS

% lower bounds
lpLB = cellfun(@(M)M.lb(bioRxnIdx)',models,'UniformOutput',false);
deltaLB = zeros(NDELTA,1);
lpLB = [[lpLB{:}]'; deltaLB];

% upper bounds
lpUB = cellfun(@(M)M.ub(bioRxnIdx)',models,'UniformOutput',false);
% minimum of allowed fold-change and cap value
deltaUB = min((epsilon-1)*minKcat,K-minKcat);
% set deltas for blacklist enzymes and unmeasured enzymes to zero
deltaUB(blackListIdx) = 0;
deltaUB(idxExclude) = 0;

lpUB = [[lpUB{:}]'; deltaUB];

% construct omega constraints
omega = sparse(2*NCOND,size(lpMatrix,2)+NCOND);
for i=1:2:2*NCOND
    condIdx = (i+1)/2;
    omega(i,[bioIdx(condIdx) size(lpMatrix,2)+condIdx]) = -[1 expVal(condIdx)];
    omega(i+1,[bioIdx(condIdx) size(lpMatrix,2)+condIdx]) = [1 -expVal(condIdx)];
end
clear condIdx

rhoMatFull = sparse(1,0);
rhoLB = zeros(0,1);
rhoUB = zeros(0,1);
rhoSumConst = sparse(0,1);
rhoSumRHS = zeros(0,1);
rhoSumSense = char;
rhoObj = zeros(0,1);

% create linear optimization problem
LP.S = [
    [lpMatrix sparse(size(lpMatrix,1),NCOND) rhoMatFull];
    [omega sparse(2*NCOND, size(rhoMatFull,2))];
    rhoSumConst];
LP.b = [
    lpRHS;
    repelem(columnVector(expVal),2,1).*repmat([-1;1],NCOND,1);
    rhoSumRHS];
LP.lb = [lpLB; zeros(NCOND,1); rhoLB];
LP.ub = [lpUB; theta.*ones(NCOND,1); rhoUB];
LP.csense = [lpSense; repmat('L',2*NCOND,1); rhoSumSense];
LP.c = [lpObj; ones(NCOND,1) / numel(models); rhoObj];

LP.osenseStr = 'min';
clear lpMatrix lpRHS lpLB lpUB lpObj omega rhoMatFull rhoLB rhoUB rhoObj ...
    rhoSumConst rhoSumRHS rhoSumSense

if any(LP.lb > LP.ub)
    fprintf('Found %d vars with lb > ub\n', sum(LP.lb > LP.ub));
end
% update reaction and metabolite names and IDs
rxnCondNrStr = strtrim(cellstr(num2str(repelem((1:NCOND)',nBiochemRxns,1))));
LP.rxns = [
    strcat(strcat(repmat(models{1}.rxns(bioRxnIdx),NCOND,1),'_cond_'), rxnCondNrStr);...
    strcat('delta_', strtok(models{1}.metNames(enzMetIdx),'['));...
    strcat({'omega_cond_'}, num2str((1:NCOND)'))
    ];

LP.rxnNames = [
    strcat(strcat(repmat(models{1}.rxnNames(bioRxnIdx),NCOND,1),{' - condition '}), rxnCondNrStr);...
    strcat('delta_', strtok(models{1}.metNames(enzMetIdx),'['));...
    strcat({'omega_cond_'}, num2str((1:NCOND)'))
    ];
if includeUM
    rxnCondNrStr = strtrim(cellstr(num2str(repelem((1:NCOND)',NDELTA,1))));
    LP.rxns = [LP.rxns;...
        strcat(repmat(strcat('rho_', strtok(models{1}.metNames(enzMetIdx),'[')),NCOND,1),...
        strcat('_cond_', rxnCondNrStr));...
        arrayfun(@(i)sprintf('delta_ptot_%d',i),1:NCOND,'un',0)'];
    LP.rxnNames = [LP.rxnNames; strcat(repmat(strcat('rho_',...
        strtok(models{1}.metNames(enzMetIdx),'[')),NCOND,1),strcat('_cond_', rxnCondNrStr));...
        arrayfun(@(i)sprintf('delta_ptot_%d',i),1:NCOND,'un',0)'];
end

metCondNrStr = strtrim(cellstr(num2str(repelem((1:NCOND)',numel(models{1}.mets),1))));
LP.mets = [strcat(strcat(repmat(models{1}.mets,NCOND,1),'_cond_'), metCondNrStr);...
    strcat(repmat({'omega_1_cond_'; 'omega_2_cond_'},NCOND,1),num2str(repelem((1:NCOND)',2,1)))] ;
LP.metNames = [strcat(strcat(repmat(models{1}.metNames,NCOND,1),{' - condition '}), metCondNrStr);...
    strcat(repmat({'omega_1_cond_'; 'omega_2_cond_'},NCOND,1),num2str(repelem((1:NCOND)',2,1)))];

if includeUM
    LP.mets = [LP.mets; strcat({'rho_sum_const_'}, num2str((1:NCOND)'))];
    LP.metNames = [LP.metNames; strcat({'rho_sum_const_'}, num2str((1:NCOND)'))];
end
clear rxnCondNrStr metCondNrStr

% solve the LP
solution = optimizeCbModel(LP);

if negCorrFlag
    
    LP_min = LP;
    
    % fix non-zero deltas and test feasibility
    deltaIdx = deltaColStartIdx:deltaColStartIdx+NDELTA-1;
    deltaVal = solution.x(deltaIdx);
    nzDeltaIdx = deltaVal~=0;
    feasTol = getCobraSolverParams('LP', 'feasTol');
    LP_min.lb(deltaIdx(nzDeltaIdx)) = deltaVal(nzDeltaIdx) - feasTol;
    LP_min.ub(deltaIdx(nzDeltaIdx)) = deltaVal(nzDeltaIdx) + feasTol;
    
    % fix uptake fluxes
    exc_idx = startsWith(LP.rxns, 'EX_');
    if sum(exc_idx) == 0
        warning('Step 2: No exchange reactions found that could be fixed.')
    end
    LP_min.lb(exc_idx) = solution.x(exc_idx) - 1e-6;
    LP_min.ub(exc_idx) = solution.x(exc_idx) + 1e-6;
    
    % set upper bound for sum of relative errors
    omegaIdx = startsWith(LP_min.rxns, 'omega_cond_');
    omegaSumRow = sparse(1, size(LP_min,2));
    omegaSumRow(omegaIdx) = 1;
    omegaSumRHS = sum(solution.x(omegaIdx)) + feasTol;
    omegaSumSense = 'L';
    LP_min.S = [LP_min.S; omegaSumRow];
    LP_min.b = [LP_min.b; omegaSumRHS];
    LP_min.csense = [LP_min.csense; omegaSumSense];
    LP_min.mets = [LP_min.mets; 'sum_omega_constraint'];
    LP_min.metNames = [LP_min.metNames; 'sum_omega_constraint'];
    
    % allow only for negative deltas if not changed in first step
    deltaLB = max((1/epsilon-1)*minKcat,min(minKcat)-minKcat);
    deltaLB(blackListIdx) = 0;
    deltaLB(idxExclude) = 0;
    
    % also exclude deltas for reactions that are inactive across all
    % conditions
    idxAllZeroFlux = false(NDELTA, 1);
    for i = 1:NDELTA
        % find reactions associated to current delta
        protID = erase(LP.rxns(deltaIdx(i)), 'delta_');
        rxnIDs = models{1}.rxns(models{1}.S(findMetIDs(models{1},protID),:)<0);
        rxnIdx = arrayfun(@(i)i+nBiochemRxns*(0:NCOND-1),...
            findRxnIDs(models{1}, rxnIDs), 'UniformOutput', false);
        rxnIdx = [rxnIdx{:}];
        idxAllZeroFlux(i) = ~any(solution.x(rxnIdx));
    end
    deltaLB(idxAllZeroFlux) = 0;
    
    LP_min.lb(deltaIdx(~nzDeltaIdx)) = deltaLB(~nzDeltaIdx);
    LP_min.ub(deltaIdx(~nzDeltaIdx)) = zeros(sum(~nzDeltaIdx),1);
    
    % also update delta lower bounds in LP that will be returned by the
    % function
    LP.lb(deltaIdx(~nzDeltaIdx)) = deltaLB(~nzDeltaIdx);
    
    % minimize the sum of allowed deltas
    LP_min.c(deltaIdx) = 0;
    LP_min.c(deltaIdx(~nzDeltaIdx)) = lambda / (sum(any(E,2)) - sum(blackListIdx) - sum(nzDeltaIdx));
    
    solution_min = optimizeCbModel(LP_min);
    
    solution = solution_min;
end

if solution.stat == 1
    % get delta values
    deltaVal = solution.x(deltaColStartIdx:deltaColStartIdx+NDELTA-1);
    % find indices of enzymes with corrected kcats
    enzMetChangedIdx = find(deltaVal~=0);
    % reduce delta values to non-zero values
    deltaVal(deltaVal==0) = [];
    enzMetNames = regexprep(LP.metNames(enzMetIdx(enzMetChangedIdx)),' - condition \d+','');
    changeTab = cell2table(repmat({'',0,0,0,0,0,0,true},numel(enzMetChangedIdx),1),'VariableNames',...
        {'ENZYME_MET_NAME','KCAT_ORIG [s^-1]','KCAT_UPDATED [s^-1]',...
        'FOLD_CHANGE','DELTA','DELTA_LB','DELTA_UB', 'MEASURED'});
    for i=1:numel(enzMetChangedIdx)
        % calculate updated kcat value
        updtKcat = minKcat(enzMetChangedIdx(i))+deltaVal(i);
        % update result table
        changeTab.(1)(i) = enzMetNames(i);
        changeTab.(2)(i) = minKcat(enzMetChangedIdx(i))/3600;
        changeTab.(3)(i) = updtKcat/3600;
        changeTab.(4)(i) = updtKcat/minKcat(enzMetChangedIdx(i));
        changeTab.(5)(i) = deltaVal(i)/3600;
        changeTab.(6)(i) = deltaLB(enzMetChangedIdx(i))/3600;
        changeTab.(7)(i) = deltaUB(enzMetChangedIdx(i))/3600;
        changeTab.(8)(i) = true;
    end
    
    if includeUM
        isMeasured = ~any(rhoMatSingle);
        changeTab.(7) = full(isMeasured(enzMetChangedIdx))';
    end
    % update kcats in models
    for i=1:numel(models)
        tmpModel = models{i};
        EC       = tmpModel.enzymeConstraints;
        kcatVec  = EC.kcat(:)*3600;
        for j=1:numel(enzMetChangedIdx)
            rowIdx = enzMetIdx(enzMetChangedIdx(j));
            enzID  = tmpModel.mets{rowIdx};
            protID = enzID(numel(enzMetPfx)+1:end);
            [~, enzCol] = ismember(protID, EC.enzymes);
            rxnRows = find(EC.rxnEnzMat(:, enzCol) ~= 0);
            origKcats = kcatVec(rxnRows);
            idxMin = origKcats==minKcat(enzMetChangedIdx(j));
            updtKcat = (minKcat(enzMetChangedIdx(j)) + deltaVal(j))/3600;
            tmpModel.enzymeConstraints.kcat(rxnRows(idxMin)) = updtKcat;
        end
        tmpModel  = UpdateSmatrix(tmpModel);
        models{i} = tmpModel;
    end
    
    % get relative error
    relError = solution.x(contains(LP.rxns,'omega_cond_'));
    
else
    changeTab = table;
    relError = zeros(sum(contains(LP.rxns,'omega_cond_')),1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = parseInput(arguments)
    % set default values
    ENZ_MET_PFX_DEFAULT = 'prot_';
    ENZ_RXN_pfx_DEFAULT = 'usage_prot_';
    EPSILON_DEFAULT = 1e5;
    THETA_DEFAULT = 0.5;
    LAMBDA_DEFAULT = 1e-7;
    K_DEFAULT = 57500000; % Pyrococcus furiosus; 5.3.1.1; D-glyceraldehyde 3-phosphate
    INCL_UM_DEFAULT = false;
    F_N_DEFAULT = 0.5;
    SIGMA_DEFAULT = 0.5;
    NEG_CORR_F_DEFAULT = false;
    
    % validation functions
    validateModel = @(M) iscell(M)||isstruct(M);
    validScalarDouble = @(v)~ischar(v)&isscalar(v);
    
    p = inputParser;
    p.FunctionName = 'PRESTO';
    
    addRequired(p,'models',validateModel)
    addRequired(p,'expVal',@isnumeric)
    addRequired(p,'E',@isnumeric)
    addParameter(p,'enzMetPfx',ENZ_MET_PFX_DEFAULT,@ischar)
    addParameter(p,'enzRxnPfx',ENZ_RXN_pfx_DEFAULT,@ischar)
    addParameter(p,'epsilon',EPSILON_DEFAULT,validScalarDouble)
    addParameter(p,'theta',THETA_DEFAULT,@iscolumn)
    addParameter(p,'lambda',LAMBDA_DEFAULT,validScalarDouble)
    addParameter(p,'enzBlackList',{''},@iscellstr)
    addParameter(p,'K',K_DEFAULT,validScalarDouble)
    addParameter(p,'includeUM',INCL_UM_DEFAULT,@islogical)
    addParameter(p,'f_n',F_N_DEFAULT,@isnumeric)
    addParameter(p,'sigma',SIGMA_DEFAULT,validScalarDouble)
    addParameter(p,'pCorrFactor',1,@isnumeric)
    addParameter(p,'negCorrFlag',NEG_CORR_F_DEFAULT,@islogical)
    
    parse(p,arguments{:})
end

function [enzMetIdx,enzRxnIdx,poolMetIdx,poolRxnIdx] = findEnzIdx(model,enzMetPfx,enzRxnPfx)
    % Find indices of enzyme metabolites and associated enzyme usage
    
    enzMetIdx = find(contains(model.mets,enzMetPfx));
        poolMetIdx = false(numel(model.mets),1);
    if isfield(model,'mets')
        poolMetIdx = poolMetIdx | strcmp(model.mets, 'prot_pool');
        poolMetIdx = find(poolMetIdx);
    end

    % Remove the line of the protein pool from the enzyme metabolite collection
    enzMetIdx = setdiff(enzMetIdx, poolMetIdx);

    if isempty(enzMetIdx)
        error('No enzyme metabolites found in model with given prefix')
    end
    
    enzRxnIdx = find(contains(model.rxns,enzRxnPfx));
    poolRxnIdx = find(contains(model.rxns,'prot_pool_exchange'));
    if sum(enzRxnIdx) == 0
        error('No enzyme usage pseudoreactions found in model with given prefix')
    end
end

function bioIdx = findBioIdx(models,bioRxnIdx, nBiochemRxns)
    % find biomass reaction indices from model objectives and return a
    % double array that contains the positions within the complete linear
    % optimization problem
    bioIdx = nan(numel(models),1);
    for m=1:numel(models)
        j_global = find(models{m}.c);
        if numel(j_global) ~= 1
            error('Model %d must have exactly one objective (biomass) reaction.', m);
        end
        j_local = find(bioRxnIdx == j_global, 1);
        bioIdx(m) = (m-1)*nBiochemRxns + j_local;
    end
    
end

function ProteinpoolIdx = findproteinpoolIdx(models,poolIdx, nBiochemRxns)
    ProteinpoolIdx = nan(numel(models),1);
    for m=1:numel(models)
        ProteinpoolIdx(m) = (m-1)*nBiochemRxns+poolIdx;
    end
end

function [origKcat,minKcat,nUniqPerEnzyme] = getKcatFromEcModel(model,enzMetPfx,enzRxnPfx)
%GETKCATFROMECMODEL  Extract per-enzyme kcat statistics from an ecModel.

    % Find indices of enzyme pseudo-metabolites in model.mets
    % (enzRxnPfx is not used directly here, but required by findEnzIdx)
    [enzMetIdx,~,~] = findEnzIdx(model,enzMetPfx,enzRxnPfx);
    nEnz = numel(enzMetIdx);

    % Short handle to enzymeConstraints
    EC = model.enzymeConstraints;
    % Flatten kcat to a column vector (one kcat per EC reaction) 
    kcatVec = EC.kcat(:)*3600; % change to h-1

    % Preallocate outputs
    minKcat        = nan(nEnz, 1);      % per-enzyme minimum kcat
    nUniqPerEnzyme = zeros(nEnz, 1);    % per-enzyme #unique kcat
    origKcat_cell  = cell(nEnz, 1);     % store kcat list for each enzyme

    % Loop over each enzyme pseudo-metabolite
    for e = 1:nEnz
        % Get metabolite ID, e.g. 'prot_P0A858_c'
        metID = model.mets{enzMetIdx(e)};
        % Strip the prefix (e.g. 'prot_') to obtain the enzyme ID
        % as stored in EC.enzymes. Here we assume no extra suffix
        % beyond the prefix (adjust if you have compartment suffixes).
        protID = metID(numel(enzMetPfx)+1:end);

        % Map this enzyme ID to a column in EC.enzymes
        [hit, enzCol] = ismember(protID, EC.enzymes);
        if ~hit
            % If no mapping is found, skip and warn
            warning('getKcatFromEcModel:NoMapping', ...
                    'Could not map enzyme metabolite "%s" (protID = "%s") to enzymeConstraints.enzymes.', ...
                    metID, protID);
            continue;
        end

        % Find all EC reactions that use this enzyme (non-zero entries
        % in the corresponding column of rxnEnzMat)
        rxnRows = find(EC.rxnEnzMat(:, enzCol) ~= 0);
        if isempty(rxnRows)
            % Enzyme is present but not connected to any EC reaction
            continue;
        end

        % Collect kcat values for those reactions
        kVals = kcatVec(rxnRows);
        % Keep only finite, positive kcat values
        kVals = kVals(isfinite(kVals) & kVals > 0);
        if isempty(kVals)
            continue;
        end

        % Store raw kcat list for this enzyme
        origKcat_cell{e}  = kVals(:);
        % Count unique kcat values used by this enzyme
        nUniqPerEnzyme(e) = numel(unique(kVals));
        % Minimum kcat for this enzyme (used by PRESTO-style formulations)
        minKcat(e)        = min(kVals);
    end

    % Concatenate all kcat values into a single vector
    if any(~cellfun(@isempty, origKcat_cell))
        origKcat = vertcat(origKcat_cell{:});
    else
        origKcat = [];
    end
end
