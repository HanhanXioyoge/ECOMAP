%% GAUKS growth-experiment validation using BayesianGrowthRates.tsv
% Dynamically discovers model pairs; no fixed model-pair counts are assumed.
% Groups: All, AboveUmin (> Umin), BelowUmin (<= Umin).
% Metrics: direct biomass RMSE, pooled exchange RMSE, and pooled total RMSE.
% The condition-specific substrate exchange is fixed to experimental Uptake and
% excluded from exchange scoring. Total RMSE mixes units; report components too.
% Unknown/limited oxygen conditions are retained in condition output but excluded
% from primary RMSE summaries until a numeric oxygen bound is specified.
% Umin=1000 denotes an inactive/unset substrate threshold: exclude that
% model-pair/substrate combination from all scores, retaining audit rows.
%
% Output workbook sheets:
%   ConditionPredictions (long format: biomass and exchange rows)
%   PairSummary, SpeciesSummary, DataAudit
%% Load models
eciML1515_Location = fullfile(findECOMAProot,'projects','eciML1515');
eciML1515_models_Location = fullfile(eciML1515_Location,'models');
eciML1515_kcatRepo_Location = fullfile(eciML1515_Location,'calibration');
load(fullfile(eciML1515_models_Location,'eciML1515_models.mat'))

ecYeast_Location = fullfile(findECOMAProot,'projects','ecYeast');
ecYeast_models_Location = fullfile(ecYeast_Location,'models');
ecYeast_kcatRepo_Location = fullfile(ecYeast_Location,'calibration');
load(fullfile(ecYeast_models_Location,'ecYeast_models.mat'))

eciCW773_Location = fullfile(findECOMAProot,'projects','eciCW773');
eciCW773_models_Location = fullfile(eciCW773_Location,'models');
eciCW773_kcatRepo_Location = fullfile(eciCW773_Location,'calibration');
load(fullfile(eciCW773_models_Location,'eciCW773_models.mat'))
%% Species/model configuration
speciesNames = ["eciML1515","ecYeast","eciCW773"];
orgNames = {'Escherichia coli','Saccharomyces cerevisiae', ...
            'Corynebacterium glutamicum'};
modelSets = {eciML1515_models,ecYeast_models,eciCW773_models};
modelLocations = {eciML1515_kcatRepo_Location, ...
                  ecYeast_kcatRepo_Location,eciCW773_kcatRepo_Location};
ecGEM_bioRxns = {'biomass','r_2111','biomass'};
oxygenRxns = {'EX_o2_e','r_1992','EX_o2_e'};
%% Read experimental data
nSpecies = numel(speciesNames);
expSets = cell(1,nSpecies);
for s = 1:nSpecies
    expSets{s} = readtable( ...
        fullfile(modelLocations{s},'data','BayesianGrowthRates.tsv'), ...
        'FileType','text','Delimiter','\t','VariableNamingRule','preserve');
end
%% Parameters and result buffers
inactiveUmin = 1000;
inactiveUminTol = 1e-6; % Absolute tolerance for the sentinel value only.
nBootPair = 2000;
nBootSpecies = 5000;
rng(1);
RCondition = cell(0,18);
RFlux = cell(0,18);
RAudit = cell(0,19);
%% Predict each discovered base/GAUKS model pair
for s = 1:nSpecies
    M = modelSets{s};
    Traw = expSets{s};
    bioRxn = ecGEM_bioRxns{s};
    assert(ismember(bioRxn,Traw.Properties.VariableNames), ...
        'GAUKS:MissingBiomassColumn', ...
        'Experimental biomass column %s is missing for %s.', ...
        bioRxn,char(speciesNames(s)));
    targetRaw = toNumeric(Traw.(bioRxn));
    uptakeRaw = toNumeric(Traw.Uptake);
    substrateRaw = string(Traw.Substrate);
    oxRaw = normalizeOxygen(Traw.OxAvail);
    % Biomass may be missing while exchange observations remain usable.
    valid = isfinite(uptakeRaw) & uptakeRaw < 0 & ...
            ~ismissing(substrateRaw) & substrateRaw ~= "" & ...
            substrateRaw ~= "NaN";
    T = Traw(valid,:);
    target = targetRaw(valid);
    uptakeAll = uptakeRaw(valid);
    conditionID = find(valid);
    oxValid = oxRaw(valid);
    nRaw = height(Traw);
    nValid = height(T);
    supportedOxygen = oxValid == "aerobic" | oxValid == "anaerobic";
    nUnsupportedOxygen = sum(~supportedOxygen);
    allCarbon = unique(string(T.Substrate),'stable');
    allCarbon = allCarbon(~ismissing(allCarbon) & ...
                          allCarbon ~= "" & allCarbon ~= "NaN");
    gauksFields = fieldnames(M);
    gauksFields = gauksFields(endsWith(gauksFields,'_GAUKS'));
    for p = 1:numel(gauksFields)
        gauksName = gauksFields{p};
        baseName = erase(gauksName,'_GAUKS');
        if ~isfield(M,baseName)
            RAudit(end+1,:) = {speciesNames(s),string(baseName), ...
                string(gauksName),"MissingBaseModel",nRaw,nValid, ...
                nUnsupportedOxygen,NaN,NaN,NaN,0,0,0,0,0,0,0,0,0};
            continue
        end
        mB = M.(baseName);
        mG = M.(gauksName);
        ibB = findRxnIDs(mB,bioRxn);
        ibG = findRxnIDs(mG,bioRxn);
        if ibB == 0 || ibG == 0
            RAudit(end+1,:) = {speciesNames(s),string(baseName), ...
                string(gauksName),"MissingBiomassReaction",nRaw,nValid, ...
                nUnsupportedOxygen,NaN,NaN,NaN,0,0,0,0,0,0,0,0,0};
            continue
        end
        mB.c(:) = 0; mG.c(:) = 0;
        mB.c(ibB) = 1; mG.c(ibG) = 1;
        idsB = findRxnIDs(mB,cellstr(allCarbon));
        idsG = findRxnIDs(mG,cellstr(allCarbon));
        params = extractSluiceParams(mG,cellstr(allCarbon));
        % Data columns that are actual reaction IDs and have numeric observations.
        varNames = Traw.Properties.VariableNames;
        % Biomass IDs are species-specific (e.g. r_2111 for ecYeast).
        % Use the same configured ID for experimental data and model objective.
        metadataNames = ["Substrate","Uptake","OxAvail","Media", ...
                         string(ecGEM_bioRxns{s})];
        fluxRxns = {};
        for v = 1:numel(varNames)
            rxn = varNames{v};
            if any(string(rxn) == metadataNames)
                continue
            end
            expValues = toNumeric(Traw.(rxn));
            if ~any(isfinite(expValues(valid)))
                continue
            end
            if findRxnIDs(mB,rxn) > 0 || findRxnIDs(mG,rxn) > 0
                fluxRxns{end+1,1} = rxn;
            end
        end
        nFluxRxns = numel(fluxRxns);
        fluxExpByRxn = cell(nFluxRxns,1);
        fluxIdsB = zeros(nFluxRxns,1);
        fluxIdsG = zeros(nFluxRxns,1);
        for r = 1:nFluxRxns
            fluxExpByRxn{r} = toNumeric(T.(fluxRxns{r}));
            fluxIdsB(r) = findRxnIDs(mB,fluxRxns{r});
            fluxIdsG(r) = findRxnIDs(mG,fluxRxns{r});
        end

        nBaseSolved = 0; nGAUKSSolved = 0;
        nBaseFailed = 0; nGAUKSFailed = 0;
        nPairedBio = 0; nPairedExchangeObs = 0;
        groupLabels = strings(height(T),1);

        %% Evaluate conditions for this model pair
        for i = 1:height(T)

            ex = string(T.Substrate(i));
            ox = oxValid(i);
            media = string(T.Media(i));
            muExp = target(i);
            uptakeExp = uptakeAll(i);
            uptakeMagnitude = abs(uptakeExp);

            k = find(allCarbon == ex,1);
            umin = NaN;
            if ~isempty(k) && isfield(params,'umin') && k <= numel(params.umin)
                umin = double(params.umin(k));
            end

            isInactiveUmin = isfinite(umin) && ...
                abs(umin-inactiveUmin) <= inactiveUminTol;
            if isInactiveUmin
                uminGroup = "ExcludedInactiveUmin";
            elseif ~isfinite(umin) || umin < 0
                uminGroup = "UnknownUmin";
            elseif uptakeMagnitude > umin
                uminGroup = "AboveUmin";
            else
                uminGroup = "BelowUmin";
            end
            groupLabels(i) = uminGroup;

            muB = NaN; muG = NaN;
            vB = []; vG = [];
            statusB = "NotSolved"; statusG = "NotSolved";

            jB = findRxnIDs(mB,char(ex));
            jG = findRxnIDs(mG,char(ex));
            iBasal = findRxnIDs(mG,char(ex+"_basal"));
            iExtended = findRxnIDs(mG,char(ex+"_extended"));

            % Do not treat unspecified "limited" oxygen as aerobic.
            if isInactiveUmin
                % Exclude both models on the same substrate, not the whole
                % model pair: other substrates may have active thresholds.
                statusB = "Excluded_InactiveUmin";
                statusG = "Excluded_InactiveUmin";
            elseif ox ~= "aerobic" && ox ~= "anaerobic"
                statusB = "Excluded_OxygenCondition";
                statusG = "Excluded_OxygenCondition";
            else
                %% Base model
                if jB == 0
                    statusB = "MissingSubstrateRxn";
                else
                    b = mB;
                    idxB = idsB(idsB > 0);
                    b.lb(idxB) = 0;
                    b.ub(idxB) = max(0,b.ub(idxB));
                    b.lb(jB) = uptakeExp;
                    b.ub(jB) = uptakeExp;

                    if ox == "anaerobic"
                        b = anaerobicModel(b,orgNames{s},oxygenRxns{s});
                    end

                    sb = solveLP(b,'max');
                    if ~isempty(sb.x) && sb.stat == 1 && isfinite(sb.x(ibB))
                        vB = sb.x;
                        muB = vB(ibB);
                        statusB = "Solved";
                    else
                        statusB = "SolveFailed";
                    end
                end

                %% GAUKS model
                if jG == 0 || iBasal == 0 || iExtended == 0
                    statusG = "MissingGAUKSRxn";
                else
                    g = mG;
                    idxG = idsG(idsG > 0);
                    g.lb(idxG) = 0;
                    g.ub(idxG) = max(0,g.ub(idxG));
                    g.lb(jG) = uptakeExp;
                    g.ub(jG) = uptakeExp;

                    if ox == "anaerobic"
                        g = anaerobicModel(g,orgNames{s},oxygenRxns{s});
                    end

                    sg = solveLP(g,'max');
                    if ~isempty(sg.x) && sg.stat == 1 && isfinite(sg.x(ibG))
                        vG = sg.x;
                        muG = vG(ibG);
                        statusG = "Solved";
                    else
                        statusG = "SolveFailed";
                    end
                end
            end

            nBaseSolved = nBaseSolved + double(statusB == "Solved");
            nGAUKSSolved = nGAUKSSolved + double(statusG == "Solved");
            nBaseFailed = nBaseFailed + double( ...
                statusB ~= "Solved" && ~startsWith(statusB,"Excluded_"));
            nGAUKSFailed = nGAUKSFailed + double( ...
                statusG ~= "Solved" && ~startsWith(statusG,"Excluded_"));
            nPairedBio = nPairedBio + double( ...
                isfinite(muExp) && isfinite(muB) && isfinite(muG));

            RCondition(end+1,:) = { ...
                speciesNames(s),string(baseName),string(gauksName), ...
                conditionID(i),ex,ox,media,uptakeExp,uptakeMagnitude, ...
                muExp,umin,uminGroup,muB,muG,muB-muExp,muG-muExp, ...
                statusB,statusG};

            %% Record other measured exchange reactions in long format
            for r = 1:nFluxRxns
                rxn = string(fluxRxns{r});

                % This condition's substrate uptake is fixed to experiment.
                if rxn == ex
                    continue
                end

                expFlux = fluxExpByRxn{r}(i);
                predFB = NaN; predFG = NaN;
                if ~isempty(vB) && fluxIdsB(r) > 0
                    predFB = vB(fluxIdsB(r));
                end
                if ~isempty(vG) && fluxIdsG(r) > 0
                    predFG = vG(fluxIdsG(r));
                end

                if isInactiveUmin
                    fluxStatusB = "Excluded_InactiveUmin";
                elseif fluxIdsB(r) == 0
                    fluxStatusB = "MissingReactionInModel";
                elseif isempty(vB)
                    fluxStatusB = "ModelNotSolved";
                elseif ~isfinite(expFlux)
                    fluxStatusB = "NoExperimentalObservation";
                else
                    fluxStatusB = "Predicted";
                end

                if isInactiveUmin
                    fluxStatusG = "Excluded_InactiveUmin";
                elseif fluxIdsG(r) == 0
                    fluxStatusG = "MissingReactionInModel";
                elseif isempty(vG)
                    fluxStatusG = "ModelNotSolved";
                elseif ~isfinite(expFlux)
                    fluxStatusG = "NoExperimentalObservation";
                else
                    fluxStatusG = "Predicted";
                end

                errFB = NaN; errFG = NaN;
                if isfinite(expFlux) && isfinite(predFB)
                    errFB = predFB-expFlux;
                end
                if isfinite(expFlux) && isfinite(predFG)
                    errFG = predFG-expFlux;
                end

                nPairedExchangeObs = nPairedExchangeObs + double( ...
                    isfinite(expFlux) && isfinite(predFB) && isfinite(predFG));

                RFlux(end+1,:) = { ...
                    speciesNames(s),string(baseName),string(gauksName), ...
                    conditionID(i),ex,ox,media,uptakeExp,umin,uminGroup, ...
                    rxn,expFlux,predFB,predFG,errFB,errFG, ...
                    fluxStatusB,fluxStatusG};
            end
        end

        %% Audit this model pair
        RAudit(end+1,:) = { ...
            speciesNames(s),string(baseName),string(gauksName),"Analyzed", ...
            nRaw,nValid,nUnsupportedOxygen, ...
            sum(groupLabels == "AboveUmin"), ...
            sum(groupLabels == "BelowUmin"), ...
            sum(groupLabels == "UnknownUmin"), ...
            nBaseSolved,nGAUKSSolved,nBaseFailed,nGAUKSFailed,nPairedBio, ...
            nFluxRxns,sum(fluxIdsB == 0),sum(fluxIdsG == 0), ...
            nPairedExchangeObs};

    end
end

%% Convert condition-level results to typed tables
conditionResults = cell2table(RCondition,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','ConditionID','Substrate', ...
    'OxAvail','Media','UptakeExperimental','UptakeMagnitude', ...
    'GrowthExperimental','Umin','UminGroup','GrowthBefore','GrowthAfter', ...
    'ErrorBefore','ErrorAfter','BaseStatus','GAUKSStatus'});

conditionResults = castTableColumns(conditionResults, ...
    {'Species','BaseModel','GAUKSModel','Substrate','OxAvail','Media', ...
     'UminGroup','BaseStatus','GAUKSStatus'}, ...
    {'ConditionID','UptakeExperimental','UptakeMagnitude', ...
     'GrowthExperimental','Umin','GrowthBefore','GrowthAfter', ...
     'ErrorBefore','ErrorAfter'});

fluxResults = cell2table(RFlux,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','ConditionID','Substrate','OxAvail', ...
    'Media','UptakeExperimental','Umin','UminGroup','Reaction', ...
    'FluxExperimental','FluxBefore','FluxAfter','ErrorBefore','ErrorAfter', ...
    'StatusBefore','StatusAfter'});

fluxResults = castTableColumns(fluxResults, ...
    {'Species','BaseModel','GAUKSModel','Substrate','OxAvail','Media', ...
     'UminGroup','Reaction','StatusBefore','StatusAfter'}, ...
    {'ConditionID','UptakeExperimental','Umin','FluxExperimental', ...
     'FluxBefore','FluxAfter','ErrorBefore','ErrorAfter'});

dataAudit = cell2table(RAudit,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','PairStatus','NRawConditions', ...
    'NValidInputs','NUnsupportedOxygen','NAboveUmin','NBelowUmin', ...
    'NUnknownUmin','NBaseSolved','NGAUKSSolved','NBaseFailed', ...
    'NGAUKSFailed','NPairedBio','NFluxReactions', ...
    'NFluxMissingInBase','NFluxMissingInGAUKS','NPairedExchangeObs'});

dataAudit = castTableColumns(dataAudit, ...
    {'Species','BaseModel','GAUKSModel','PairStatus'}, ...
    {'NRawConditions','NValidInputs','NUnsupportedOxygen','NAboveUmin', ...
     'NBelowUmin','NUnknownUmin','NBaseSolved','NGAUKSSolved', ...
     'NBaseFailed','NGAUKSFailed','NPairedBio','NFluxReactions', ...
     'NFluxMissingInBase','NFluxMissingInGAUKS','NPairedExchangeObs'});

% Keep inactive-threshold exclusions distinct from solver failures.
dataAudit.NExcludedInactiveUmin = zeros(height(dataAudit),1);
dataAudit.NEligibleConditions = zeros(height(dataAudit),1);
for p = 1:height(dataAudit)
    pairRows = conditionResults.Species == dataAudit.Species(p) & ...
        conditionResults.BaseModel == dataAudit.BaseModel(p) & ...
        conditionResults.GAUKSModel == dataAudit.GAUKSModel(p);
    dataAudit.NExcludedInactiveUmin(p) = sum(pairRows & ...
        conditionResults.UminGroup == "ExcludedInactiveUmin");
    dataAudit.NEligibleConditions(p) = sum(pairRows & ...
        conditionResults.UminGroup ~= "ExcludedInactiveUmin" & ...
        (conditionResults.OxAvail == "aerobic" | ...
         conditionResults.OxAvail == "anaerobic"));
end

%% Per-model-pair RMSE summary: All / AboveUmin / BelowUmin
pairKeys = unique(conditionResults(:, ...
    {'Species','BaseModel','GAUKSModel'}),'rows','stable');
groupNames = ["All","AboveUmin","BelowUmin"];
RPairSummary = cell(0,23);

for p = 1:height(pairKeys)

    maskPairC = conditionResults.Species == pairKeys.Species(p) & ...
                conditionResults.BaseModel == pairKeys.BaseModel(p) & ...
                conditionResults.GAUKSModel == pairKeys.GAUKSModel(p);
    maskPairF = fluxResults.Species == pairKeys.Species(p) & ...
                fluxResults.BaseModel == pairKeys.BaseModel(p) & ...
                fluxResults.GAUKSModel == pairKeys.GAUKSModel(p);

    for gg = 1:numel(groupNames)
        groupName = groupNames(gg);

        maskGroupC = groupMask(conditionResults.UminGroup,groupName);
        supportedC = conditionResults.OxAvail == "aerobic" | ...
                     conditionResults.OxAvail == "anaerobic";
        C = conditionResults(maskPairC & maskGroupC & supportedC,:);

        maskGroupF = groupMask(fluxResults.UminGroup,groupName);
        supportedF = fluxResults.OxAvail == "aerobic" | ...
                     fluxResults.OxAvail == "anaerobic";
        F = fluxResults(maskPairF & maskGroupF & supportedF,:);

        nConditions = height(C);
        pairedBio = isfinite(C.GrowthExperimental) & ...
                    isfinite(C.GrowthBefore) & isfinite(C.GrowthAfter);
        nPairedBio = sum(pairedBio);
        bioResidualB = C.GrowthBefore(pairedBio)-C.GrowthExperimental(pairedBio);
        bioResidualG = C.GrowthAfter(pairedBio)-C.GrowthExperimental(pairedBio);

        rmseBioB = NaN; rmseBioG = NaN; deltaBio = NaN;
        if nPairedBio > 0
            rmseBioB = sqrt(mean(bioResidualB.^2));
            rmseBioG = sqrt(mean(bioResidualG.^2));
            deltaBio = rmseBioG-rmseBioB;
        end

        pairedFlux = isfinite(F.FluxExperimental) & ...
                     isfinite(F.FluxBefore) & isfinite(F.FluxAfter);
        nPairedExchange = sum(pairedFlux);
        fluxResidualB = F.FluxBefore(pairedFlux)-F.FluxExperimental(pairedFlux);
        fluxResidualG = F.FluxAfter(pairedFlux)-F.FluxExperimental(pairedFlux);

        rmseExchangeB = NaN; rmseExchangeG = NaN; deltaExchange = NaN;
        if nPairedExchange > 0
            rmseExchangeB = sqrt(mean(fluxResidualB.^2));
            rmseExchangeG = sqrt(mean(fluxResidualG.^2));
            deltaExchange = rmseExchangeG-rmseExchangeB;
        end

        % Total RMSE pools all paired biomass and exchange residuals directly.
        nTotalObservations = nPairedBio+nPairedExchange;
        rmseTotalB = NaN; rmseTotalG = NaN; deltaTotal = NaN;
        % A combined score requires observations from both components.
        if nPairedBio > 0 && nPairedExchange > 0
            rmseTotalB = sqrt((sum(bioResidualB.^2)+sum(fluxResidualB.^2)) / ...
                              nTotalObservations);
            rmseTotalG = sqrt((sum(bioResidualG.^2)+sum(fluxResidualG.^2)) / ...
                              nTotalObservations);
            deltaTotal = rmseTotalG-rmseTotalB;
        end

        ciBioLow = NaN; ciBioHigh = NaN;
        ciExchangeLow = NaN; ciExchangeHigh = NaN;
        ciTotalLow = NaN; ciTotalHigh = NaN;
        if nConditions >= 2
            [ciBioLow,ciBioHigh,ciExchangeLow,ciExchangeHigh, ...
             ciTotalLow,ciTotalHigh] = bootstrapPairRMSEDeltas(C,F,nBootPair);
        end

        RPairSummary(end+1,:) = { ...
            pairKeys.Species(p),pairKeys.BaseModel(p),pairKeys.GAUKSModel(p), ...
            groupName,nConditions,nPairedBio,nPairedExchange,nTotalObservations, ...
            rmseBioB,rmseBioG,deltaBio,ciBioLow,ciBioHigh, ...
            rmseExchangeB,rmseExchangeG,deltaExchange, ...
            ciExchangeLow,ciExchangeHigh,rmseTotalB,rmseTotalG, ...
            deltaTotal,ciTotalLow,ciTotalHigh};
    end
end

pairSummary = cell2table(RPairSummary,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','Group','NConditions', ...
    'NPairedBio','NPairedExchange','NTotalObservations', ...
    'BioRMSEBefore','BioRMSEAfter','DeltaBioRMSE', ...
    'BioDeltaCI_Lower','BioDeltaCI_Upper', ...
    'ExchangeRMSEBefore','ExchangeRMSEAfter','DeltaExchangeRMSE', ...
    'ExchangeDeltaCI_Lower','ExchangeDeltaCI_Upper', ...
    'TotalRMSEBefore','TotalRMSEAfter','DeltaTotalRMSE', ...
    'TotalDeltaCI_Lower','TotalDeltaCI_Upper'});

pairSummary = castTableColumns(pairSummary, ...
    {'Species','BaseModel','GAUKSModel','Group'}, ...
    {'NConditions','NPairedBio','NPairedExchange','NTotalObservations', ...
     'BioRMSEBefore','BioRMSEAfter','DeltaBioRMSE', ...
     'BioDeltaCI_Lower','BioDeltaCI_Upper','ExchangeRMSEBefore', ...
     'ExchangeRMSEAfter','DeltaExchangeRMSE','ExchangeDeltaCI_Lower', ...
     'ExchangeDeltaCI_Upper','TotalRMSEBefore','TotalRMSEAfter', ...
     'DeltaTotalRMSE','TotalDeltaCI_Lower','TotalDeltaCI_Upper'});

%% Species-level summary across discovered model pairs
% Model-pair bootstrap intervals describe configuration variability only.
% Shared data/parameters mean these are not biological-replicate intervals.
metricNames = ["BioRMSE","ExchangeRMSE","TotalRMSE"];
beforeColumns = {'BioRMSEBefore','ExchangeRMSEBefore','TotalRMSEBefore'};
afterColumns = {'BioRMSEAfter','ExchangeRMSEAfter','TotalRMSEAfter'};
deltaColumns = {'DeltaBioRMSE','DeltaExchangeRMSE','DeltaTotalRMSE'};
RSpeciesSummary = cell(0,14);

for s = 1:nSpecies
    for gg = 1:numel(groupNames)
        groupName = groupNames(gg);

        for mm = 1:numel(metricNames)
            mask = pairSummary.Species == speciesNames(s) & ...
                   pairSummary.Group == groupName;

            beforeAll = pairSummary.(beforeColumns{mm})(mask);
            afterAll = pairSummary.(afterColumns{mm})(mask);
            deltaAll = pairSummary.(deltaColumns{mm})(mask);
            validPair = isfinite(beforeAll) & isfinite(afterAll) & ...
                        isfinite(deltaAll);

            beforeValues = beforeAll(validPair);
            afterValues = afterAll(validPair);
            deltaValues = deltaAll(validPair);

            nPairsTotal = sum(mask);
            nPairsWithResult = numel(deltaValues);
            meanBefore = NaN; meanAfter = NaN; meanDelta = NaN;
            medianDelta = NaN; q1Delta = NaN; q3Delta = NaN;
            improvedFraction = NaN; ciLow = NaN; ciHigh = NaN;

            if nPairsWithResult > 0
                meanBefore = mean(beforeValues);
                meanAfter = mean(afterValues);
                meanDelta = mean(deltaValues);
                medianDelta = median(deltaValues);
                q1Delta = simplePercentile(deltaValues,25);
                q3Delta = simplePercentile(deltaValues,75);
                improvedFraction = mean(deltaValues < 0);
                [ciLow,ciHigh] = bootstrapMeanCI(deltaValues,nBootSpecies);
            end

            RSpeciesSummary(end+1,:) = { ...
                speciesNames(s),groupName,metricNames(mm), ...
                nPairsTotal,nPairsWithResult,meanBefore,meanAfter, ...
                meanDelta,medianDelta,q1Delta,q3Delta, ...
                improvedFraction,ciLow,ciHigh};
        end
    end
end

speciesSummary = cell2table(RSpeciesSummary,'VariableNames',{ ...
    'Species','Group','Metric','NModelPairsTotal','NModelPairsWithResult', ...
    'MeanRMSEBefore','MeanRMSEAfter','MeanDelta','MedianDelta', ...
    'Q1Delta','Q3Delta','ImprovedFraction', ...
    'ConfigBootstrapCI_Lower','ConfigBootstrapCI_Upper'});

speciesSummary = castTableColumns(speciesSummary, ...
    {'Species','Group','Metric'}, ...
    {'NModelPairsTotal','NModelPairsWithResult','MeanRMSEBefore', ...
     'MeanRMSEAfter','MeanDelta','MedianDelta','Q1Delta','Q3Delta', ...
     'ImprovedFraction','ConfigBootstrapCI_Lower','ConfigBootstrapCI_Upper'});

%% Build a single long-format ConditionPredictions table
% Each condition has a biomass row plus one row for each scored exchange
% reaction. MetricType distinguishes the two kinds of observations.
nBioRows = height(conditionResults);
bioReaction = strings(nBioRows,1);
for s = 1:nSpecies
    bioReaction(conditionResults.Species == speciesNames(s)) = ...
        string(ecGEM_bioRxns{s});
end

bioOutput = table( ...
    conditionResults.Species,conditionResults.BaseModel, ...
    conditionResults.GAUKSModel,conditionResults.ConditionID, ...
    conditionResults.Substrate,conditionResults.OxAvail, ...
    conditionResults.Media,conditionResults.UptakeExperimental, ...
    conditionResults.UptakeMagnitude,conditionResults.Umin, ...
    conditionResults.UminGroup,repmat("Biomass",nBioRows,1), ...
    bioReaction,conditionResults.GrowthExperimental, ...
    conditionResults.GrowthBefore,conditionResults.GrowthAfter, ...
    conditionResults.ErrorBefore,conditionResults.ErrorAfter, ...
    conditionResults.BaseStatus,conditionResults.GAUKSStatus, ...
    'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','ConditionID','Substrate', ...
    'OxAvail','Media','UptakeExperimental','UptakeMagnitude','Umin', ...
    'UminGroup','MetricType','Reaction','ExperimentalValue', ...
    'PredictionBefore','PredictionAfter','ErrorBefore','ErrorAfter', ...
    'StatusBefore','StatusAfter'});

nFluxRows = height(fluxResults);
fluxOutput = table( ...
    fluxResults.Species,fluxResults.BaseModel,fluxResults.GAUKSModel, ...
    fluxResults.ConditionID,fluxResults.Substrate,fluxResults.OxAvail, ...
    fluxResults.Media,fluxResults.UptakeExperimental, ...
    abs(fluxResults.UptakeExperimental),fluxResults.Umin, ...
    fluxResults.UminGroup,repmat("Exchange",nFluxRows,1), ...
    fluxResults.Reaction,fluxResults.FluxExperimental, ...
    fluxResults.FluxBefore,fluxResults.FluxAfter, ...
    fluxResults.ErrorBefore,fluxResults.ErrorAfter, ...
    fluxResults.StatusBefore,fluxResults.StatusAfter, ...
    'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','ConditionID','Substrate', ...
    'OxAvail','Media','UptakeExperimental','UptakeMagnitude','Umin', ...
    'UminGroup','MetricType','Reaction','ExperimentalValue', ...
    'PredictionBefore','PredictionAfter','ErrorBefore','ErrorAfter', ...
    'StatusBefore','StatusAfter'});

conditionPredictions = [bioOutput;fluxOutput];
conditionPredictions.IncludedInPairedRMSE = ...
    isfinite(conditionPredictions.ExperimentalValue) & ...
    isfinite(conditionPredictions.PredictionBefore) & ...
    isfinite(conditionPredictions.PredictionAfter);

%% Export the four requested result sheets
outDir = fileparts(mfilename('fullpath'));
if ~exist(outDir,'dir')
    mkdir(outDir);
end

outFile = fullfile(outDir,'GAUKS_Growth_Experiment_RMSE_Analysis.xlsx');
if exist(outFile,'file')
    delete(outFile);
end

writetable(conditionPredictions,outFile,'Sheet','ConditionPredictions');
writetable(pairSummary,outFile,'Sheet','PairSummary');
writetable(speciesSummary,outFile,'Sheet','SpeciesSummary');
writetable(dataAudit,outFile,'Sheet','DataAudit');

disp("Growth experiment RMSE analysis complete.");
disp("Discovered GAUKS model pairs: " + string(height(dataAudit)));
disp("Analyzed model pairs: " + string(sum(dataAudit.PairStatus == "Analyzed")));
disp("Results saved to: " + string(outFile));

%% Local functions

function y = toNumeric(x)
    if isnumeric(x) || islogical(x)
        y = double(x(:));
    else
        y = str2double(string(x(:)));
    end
end

function ox = normalizeOxygen(x)
    ox = lower(strtrim(string(x(:))));
    for i = 1:numel(ox)
        % Check anaerobic first because "anaerobic" contains "aerobic".
        if contains(ox(i),"anaerob")
            ox(i) = "anaerobic";
        elseif contains(ox(i),"aerob")
            ox(i) = "aerobic";
        end
    end
end

function T = castTableColumns(T,stringVars,numericVars)
    for i = 1:numel(stringVars)
        v = stringVars{i};
        T.(v) = string(T.(v));
    end

    for i = 1:numel(numericVars)
        v = numericVars{i};
        x = T.(v);
        if iscell(x)
            if isempty(x)
                T.(v) = zeros(0,1);
            else
                T.(v) = cellfun(@(z) double(z),x);
            end
        else
            T.(v) = double(x);
        end
    end
end

function mask = groupMask(groupColumn,groupName)
    if groupName == "All"
        mask = groupColumn ~= "ExcludedInactiveUmin";
    else
        mask = groupColumn == groupName;
    end
end

function value = simplePercentile(x,p)
    x = sort(x(isfinite(x)));
    if isempty(x)
        value = NaN;
        return
    end
    if numel(x) == 1
        value = x(1);
        return
    end

    position = 1+(numel(x)-1)*p/100;
    lo = floor(position);
    hi = ceil(position);
    if lo == hi
        value = x(lo);
    else
        value = x(lo)+(position-lo)*(x(hi)-x(lo));
    end
end

function [ciLow,ciHigh] = percentileCI(x)
    x = x(isfinite(x));
    if isempty(x)
        ciLow = NaN;
        ciHigh = NaN;
        return
    end
    ciLow = simplePercentile(x,2.5);
    ciHigh = simplePercentile(x,97.5);
end

function [ciBioLow,ciBioHigh,ciExchangeLow,ciExchangeHigh, ...
          ciTotalLow,ciTotalHigh] = bootstrapPairRMSEDeltas(C,F,nBoot)

    ciBioLow = NaN; ciBioHigh = NaN;
    ciExchangeLow = NaN; ciExchangeHigh = NaN;
    ciTotalLow = NaN; ciTotalHigh = NaN;

    nConditions = height(C);
    if nConditions < 2
        return
    end

    bootBio = NaN(nBoot,1);
    bootExchange = NaN(nBoot,1);
    bootTotal = NaN(nBoot,1);

    % Cache per-condition sums/counts; resample whole conditions as clusters.
    bioOK = isfinite(C.GrowthExperimental) & ...
            isfinite(C.GrowthBefore) & isfinite(C.GrowthAfter);
    bioSqB = zeros(nConditions,1); bioSqG = zeros(nConditions,1);
    bioSqB(bioOK) = (C.GrowthBefore(bioOK)-C.GrowthExperimental(bioOK)).^2;
    bioSqG(bioOK) = (C.GrowthAfter(bioOK)-C.GrowthExperimental(bioOK)).^2;
    fluxOK = isfinite(F.FluxExperimental) & ...
             isfinite(F.FluxBefore) & isfinite(F.FluxAfter);
    [inGroup,conditionIndex] = ismember(F.ConditionID,C.ConditionID);
    fluxOK = fluxOK & inGroup;
    fluxCount = accumarray(conditionIndex(fluxOK),ones(sum(fluxOK),1),[nConditions,1]);
    fluxSqB = accumarray(conditionIndex(fluxOK), ...
        (F.FluxBefore(fluxOK)-F.FluxExperimental(fluxOK)).^2,[nConditions,1]);
    fluxSqG = accumarray(conditionIndex(fluxOK), ...
        (F.FluxAfter(fluxOK)-F.FluxExperimental(fluxOK)).^2,[nConditions,1]);

    for b = 1:nBoot
        idx = randi(nConditions,nConditions,1);
        nBio = sum(bioOK(idx)); nExchange = sum(fluxCount(idx));
        sbBio = sum(bioSqB(idx)); sgBio = sum(bioSqG(idx));
        sbFlux = sum(fluxSqB(idx)); sgFlux = sum(fluxSqG(idx));
        if nBio > 0
            bootBio(b) = sqrt(sgBio/nBio)-sqrt(sbBio/nBio);
        end
        if nExchange > 0
            bootExchange(b) = sqrt(sgFlux/nExchange)-sqrt(sbFlux/nExchange);
        end
        nTotal = nBio+nExchange;
        if nBio > 0 && nExchange > 0
            bootTotal(b) = sqrt((sgBio+sgFlux)/nTotal)- ...
                           sqrt((sbBio+sbFlux)/nTotal);
        end
    end

    if sum(bioOK) >= 2
        [ciBioLow,ciBioHigh] = percentileCI(bootBio);
    end
    if sum(fluxCount > 0) >= 2
        [ciExchangeLow,ciExchangeHigh] = percentileCI(bootExchange);
    end
    if sum(bioOK | fluxCount > 0) >= 2 && any(bioOK) && any(fluxCount > 0)
        [ciTotalLow,ciTotalHigh] = percentileCI(bootTotal);
    end
end

function [ciLow,ciHigh] = bootstrapMeanCI(values,nBoot)
    values = values(isfinite(values));
    if numel(values) < 2
        ciLow = NaN;
        ciHigh = NaN;
        return
    end

    n = numel(values);
    bootMeans = NaN(nBoot,1);
    for b = 1:nBoot
        idx = randi(n,n,1);
        bootMeans(b) = mean(values(idx));
    end
    [ciLow,ciHigh] = percentileCI(bootMeans);
end
