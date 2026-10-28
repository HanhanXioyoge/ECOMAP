%% GAUKS growth experiment validation
% Compare each BaseModel/GAUKSModel pair on BayesianGrowthRates.tsv.
%
% Analysis levels:
%   1) Species
%   2) Model pair
%   3) All / AboveUmin / AtOrBelowUmin
%
% Error metrics:
%   - Biomass RMSE (primary)
%   - Per-reaction exchange-flux RMSE and NRMSE
%   - Exchange macro-average NRMSE (equal weight per reaction)
%   - Combined NRMSE = (biomass NRMSE + exchange macro NRMSE)/2
%
% Important:
%   - The experimental biomass column is "biomass" for all species.
%     ecGEM_bioRxns specifies the model reaction used as the biomass objective.
%   - Each condition's substrate exchange reaction is fixed to T.Uptake.
%     That fixed substrate reaction is excluded from the "other exchange flux"
%     error for that condition.
%   - Exchange-flux errors are only computed when an experimental value exists.
%   - Umin is read from the corresponding GAUKS model for the condition's substrate.
%   - "All" includes valid conditions regardless of Umin group.
%   - Only aerobic/anaerobic conditions are included in primary summaries.
%     "Limited" oxygen conditions are retained in the condition/audit outputs,
%     but excluded because this script does not specify a numeric O2 uptake limit.
%   - Model-pair bootstrap resamples conditions; species-level bootstrap resamples
%     model pairs. Model pairs may not be statistically independent.

%% Load models
eciML1515_Location = fullfile(findECOMAProot, 'projects', 'eciML1515');
eciML1515_models_Location = fullfile(eciML1515_Location, 'models');
eciML1515_kcatRepo_Location = fullfile(eciML1515_Location, 'calibration');
load(fullfile(eciML1515_models_Location, 'eciML1515_models.mat'))

ecYeast_Location = fullfile(findECOMAProot, 'projects', 'ecYeast');
ecYeast_models_Location = fullfile(ecYeast_Location, 'models');
ecYeast_kcatRepo_Location = fullfile(ecYeast_Location, 'calibration');
load(fullfile(ecYeast_models_Location, 'ecYeast_models.mat'))

eciCW773_Location = fullfile(findECOMAProot, 'projects', 'eciCW773');
eciCW773_models_Location = fullfile(eciCW773_Location, 'models');
eciCW773_kcatRepo_Location = fullfile(eciCW773_Location, 'calibration');
load(fullfile(eciCW773_models_Location, 'eciCW773_models.mat'))
%% Define species and model settings
speciesNames = ["eciML1515","ecYeast","eciCW773"];
orgNames = {'Escherichia coli', ...
            'Saccharomyces cerevisiae', ...
            'Corynebacterium glutamicum'};
modelSets = {eciML1515_models, ...
             ecYeast_models, ...
             eciCW773_models};
modelLocations = {eciML1515_kcatRepo_Location, ...
                  ecYeast_kcatRepo_Location, ...
                  eciCW773_kcatRepo_Location};
% Model reaction IDs used as biomass objectives.
ecGEM_bioRxns = {'biomass','r_2111','biomass'};
oxygenRxns = {'EX_o2_e','r_1992','EX_o2_e'};
%% Read experimental data
expSets = cell(1,3);
for s = 1:3
    expSets{s} = readtable( ...
        fullfile(modelLocations{s}, 'data', 'BayesianGrowthRates.tsv'), ...
        'FileType', 'text', ...
        'Delimiter', '\t', ...
        'VariableNamingRule', 'preserve');
end
%% Analysis parameters
uptakeTol = 1e-6;
nBootPair = 1000;
nBootSpecies = 5000;
rng(1);
%% Result buffers
RCondition = cell(0,18);
RFlux = cell(0,19);
RAudit = cell(0,15);
bioScales = NaN(1,numel(speciesNames));
%% Predict biomass and measured exchange fluxes
for s = 1:numel(speciesNames)
    M = modelSets{s};
    Traw = expSets{s};
    bioRxn = ecGEM_bioRxns{s};
    % Experimental biomass column is consistently named "biomass".
    targetRaw = toNumeric(Traw.('biomass'));
    uptakeRaw = toNumeric(Traw.('Uptake'));
    substrateRaw = string(Traw.Substrate);
    oxRaw = normalizeOxygen(Traw.OxAvail);
    % Include all substrate names when closing alternative carbon uptake.
    allCarbon = unique(substrateRaw,'stable');
    allCarbon = allCarbon(~ismissing(allCarbon) & ...
                          allCarbon ~= "" & ...
                          allCarbon ~= "NaN");
    % Conditions require an observed biomass value and a negative uptake rate.
    valid = isfinite(targetRaw) & isfinite(uptakeRaw) & ...
            uptakeRaw < 0 & ~ismissing(substrateRaw) & ...
            substrateRaw ~= "";
    T = Traw(valid,:);
    target = targetRaw(valid);
    uptakeAll = uptakeRaw(valid);
    conditionID = find(valid);
    oxValid = oxRaw(valid);
    % Fixed species-level scale for normalized biomass RMSE.
    bioScales(s) = calculateScale(targetRaw);
    gauksFields = fieldnames(M);
    gauksFields = gauksFields(endsWith(gauksFields,'_GAUKS'));
    for p = 1:numel(gauksFields)

        gauksName = gauksFields{p};
        baseName = erase(gauksName,'_GAUKS');

        if ~isfield(M,baseName)
            continue
        end

        mB = M.(baseName);
        mG = M.(gauksName);

        % Set biomass objectives.
        ibB = findRxnIDs(mB,bioRxn);
        ibG = findRxnIDs(mG,bioRxn);

        if ibB == 0 || ibG == 0
            continue
        end

        mB.c(:) = 0;
        mG.c(:) = 0;
        mB.c(ibB) = 1;
        mG.c(ibG) = 1;

        % Carbon exchange reaction IDs.
        idsB = findRxnIDs(mB,cellstr(allCarbon));
        idsG = findRxnIDs(mG,cellstr(allCarbon));

        % Umin parameters from this GAUKS model.
        params = extractSluiceParams(mG,cellstr(allCarbon));

        % Identify measured flux columns that correspond to reactions in either
        % model. Metadata/index columns are automatically ignored.
        varNames = Traw.Properties.VariableNames;
        metadataNames = {'Substrate','Uptake','biomass','OxAvail','Media'};
        fluxRxns = {};

        for v = 1:numel(varNames)
            rxn = varNames{v};

            if any(strcmp(rxn,metadataNames))
                continue
            end

            idB = findRxnIDs(mB,rxn);
            idG = findRxnIDs(mG,rxn);

            if idB > 0 || idG > 0
                fluxRxns{end+1,1} = rxn;
            end
        end

        % Experimental fluxes, fixed normalization scales, and reaction IDs.
        nFluxRxns = numel(fluxRxns);
        fluxExpByRxn = cell(nFluxRxns,1);
        fluxScale = NaN(nFluxRxns,1);
        fluxIdsB = zeros(nFluxRxns,1);
        fluxIdsG = zeros(nFluxRxns,1);

        for r = 1:nFluxRxns
            fluxExpByRxn{r} = toNumeric(T.(fluxRxns{r}));
            fluxScale(r) = calculateScale(toNumeric(Traw.(fluxRxns{r})));
            fluxIdsB(r) = findRxnIDs(mB,fluxRxns{r});
            fluxIdsG(r) = findRxnIDs(mG,fluxRxns{r});
        end

        nBaseSolved = 0;
        nGAUKSSolved = 0;
        nPairedBio = 0;
        groupLabels = strings(height(T),1);

        %% Evaluate each experimental condition

        for i = 1:height(T)

            ex = string(T.Substrate(i));
            ox = oxValid(i);
            media = string(T.Media(i));
            muExp = target(i);
            uptakeExp = uptakeAll(i);
            uptakeMagnitude = abs(uptakeExp);

            % Use Umin from the current GAUKS model and current substrate.
            k = find(allCarbon == ex,1);
            umin = NaN;

            if ~isempty(k) && isfield(params,'umin') && ...
                    k <= numel(params.umin)
                umin = double(params.umin(k));
            end

            % Exactly two mutually exclusive Umin groups.
            if ~isfinite(umin) || umin < 0
                uminGroup = "UnknownUmin";
            elseif uptakeMagnitude > umin + uptakeTol
                uminGroup = "AboveUmin";
            else
                uminGroup = "AtOrBelowUmin";
            end

            groupLabels(i) = uminGroup;

            muB = NaN;
            muG = NaN;
            vB = [];
            vG = [];
            statusB = "NotSolved";
            statusG = "NotSolved";

            jB = findRxnIDs(mB,char(ex));
            jG = findRxnIDs(mG,char(ex));
            iBasal = findRxnIDs(mG,char(ex+"_basal"));
            iExtended = findRxnIDs(mG,char(ex+"_extended"));

            % No numeric O2 bound is supplied for "limited" conditions.
            if ox ~= "aerobic" && ox ~= "anaerobic"

                statusB = "Excluded_OxygenCondition";
                statusG = "Excluded_OxygenCondition";

            else

                %% Base model

                if jB == 0

                    statusB = "MissingSubstrateRxn";

                else

                    b = mB;

                    % Close uptake of other listed carbon substrates.
                    idxB = idsB(idsB > 0);
                    b.lb(idxB) = 0;
                    b.ub(idxB) = max(0,b.ub(idxB));

                    % Fix the tested substrate uptake to the measured value.
                    b.lb(jB) = uptakeExp;
                    b.ub(jB) = uptakeExp;

                    if ox == "anaerobic"
                        b = anaerobicModel(b,orgNames{s},oxygenRxns{s});
                    end

                    sb = solveLP(b,'max');

                    if ~isempty(sb.x) && sb.stat == 1 && ...
                            isfinite(sb.x(ibB))
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

                    % Close uptake of other listed carbon substrates.
                    idxG = idsG(idsG > 0);
                    g.lb(idxG) = 0;
                    g.ub(idxG) = max(0,g.ub(idxG));

                    % Fix the same substrate uptake as in the base model.
                    g.lb(jG) = uptakeExp;
                    g.ub(jG) = uptakeExp;

                    if ox == "anaerobic"
                        g = anaerobicModel(g,orgNames{s},oxygenRxns{s});
                    end

                    sg = solveLP(g,'max');

                    if ~isempty(sg.x) && sg.stat == 1 && ...
                            isfinite(sg.x(ibG))
                        vG = sg.x;
                        muG = vG(ibG);
                        statusG = "Solved";
                    else
                        statusG = "SolveFailed";
                    end

                end
            end

            nBaseSolved = nBaseSolved + double(isfinite(muB));
            nGAUKSSolved = nGAUKSSolved + double(isfinite(muG));
            nPairedBio = nPairedBio + ...
                double(isfinite(muB) && isfinite(muG));

            RCondition(end+1,:) = { ...
                speciesNames(s), string(baseName), string(gauksName), ...
                conditionID(i), ex, ox, media, uptakeExp, uptakeMagnitude, ...
                muExp, umin, uminGroup, muB, muG, ...
                muB-muExp, muG-muExp, statusB, statusG};

            %% Record all measured exchange reactions except the fixed substrate

            for r = 1:nFluxRxns

                rxn = string(fluxRxns{r});

                % This condition's substrate exchange is fixed to the observed
                % Uptake value; do not count it as a predicted exchange flux.
                if rxn == ex
                    continue
                end

                expFlux = fluxExpByRxn{r}(i);
                predFB = NaN;
                predFG = NaN;

                if ~isempty(vB) && fluxIdsB(r) > 0
                    predFB = vB(fluxIdsB(r));
                end

                if ~isempty(vG) && fluxIdsG(r) > 0
                    predFG = vG(fluxIdsG(r));
                end

                if fluxIdsB(r) == 0
                    fluxStatusB = "MissingReactionInModel";
                elseif isempty(vB)
                    fluxStatusB = "ModelNotSolved";
                elseif ~isfinite(expFlux)
                    fluxStatusB = "NoExperimentalObservation";
                else
                    fluxStatusB = "Predicted";
                end

                if fluxIdsG(r) == 0
                    fluxStatusG = "MissingReactionInModel";
                elseif isempty(vG)
                    fluxStatusG = "ModelNotSolved";
                elseif ~isfinite(expFlux)
                    fluxStatusG = "NoExperimentalObservation";
                else
                    fluxStatusG = "Predicted";
                end

                errFB = NaN;
                errFG = NaN;

                if isfinite(expFlux) && isfinite(predFB)
                    errFB = predFB-expFlux;
                end

                if isfinite(expFlux) && isfinite(predFG)
                    errFG = predFG-expFlux;
                end

                RFlux(end+1,:) = { ...
                    speciesNames(s), string(baseName), string(gauksName), ...
                    conditionID(i), ex, ox, media, uptakeExp, umin, ...
                    uminGroup, rxn, expFlux, predFB, predFG, ...
                    errFB, errFG, fluxScale(r), fluxStatusB, fluxStatusG};

            end
        end

        %% Audit this model pair

        nUnsupportedOxygen = sum( ...
            oxValid ~= "aerobic" & oxValid ~= "anaerobic");

        RAudit(end+1,:) = { ...
            speciesNames(s), string(baseName), string(gauksName), ...
            height(Traw), height(T), nUnsupportedOxygen, ...
            sum(groupLabels == "AboveUmin"), ...
            sum(groupLabels == "AtOrBelowUmin"), ...
            sum(groupLabels == "UnknownUmin"), ...
            nBaseSolved, nGAUKSSolved, nPairedBio, nFluxRxns, ...
            sum(fluxIdsB == 0), sum(fluxIdsG == 0)};

    end
end

%% Convert condition-level outputs to tables

conditionResults = cell2table(RCondition,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','ConditionID', ...
    'Substrate','OxAvail','Media','UptakeExperimental', ...
    'UptakeMagnitude','GrowthExperimental','Umin','UminGroup', ...
    'GrowthBefore','GrowthAfter','ErrorBefore','ErrorAfter', ...
    'BaseStatus','GAUKSStatus'});

conditionResults = normalizeCellTable(conditionResults, ...
    {'Species','BaseModel','GAUKSModel','Substrate','OxAvail','Media', ...
     'UminGroup','BaseStatus','GAUKSStatus'}, ...
    {'ConditionID','UptakeExperimental','UptakeMagnitude', ...
     'GrowthExperimental','Umin','GrowthBefore','GrowthAfter', ...
     'ErrorBefore','ErrorAfter'});

fluxResults = cell2table(RFlux,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','ConditionID','Substrate', ...
    'OxAvail','Media','UptakeExperimental','Umin','UminGroup', ...
    'Reaction','FluxExperimental','FluxBefore','FluxAfter', ...
    'ErrorBefore','ErrorAfter','Scale','StatusBefore','StatusAfter'});

fluxResults = normalizeCellTable(fluxResults, ...
    {'Species','BaseModel','GAUKSModel','Substrate','OxAvail','Media', ...
     'UminGroup','Reaction','StatusBefore','StatusAfter'}, ...
    {'ConditionID','UptakeExperimental','Umin','FluxExperimental', ...
     'FluxBefore','FluxAfter','ErrorBefore','ErrorAfter','Scale'});

dataAudit = cell2table(RAudit,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','NRawConditions','NValidInputs', ...
    'NUnsupportedOxygen','NAboveUmin','NAtOrBelowUmin','NUnknownUmin', ...
    'NBaseSolved','NGAUKSSolved','NPairedBio','NFluxReactions', ...
    'NFluxMissingInBase','NFluxMissingInGAUKS'});

dataAudit = normalizeCellTable(dataAudit, ...
    {'Species','BaseModel','GAUKSModel'}, ...
    {'NRawConditions','NValidInputs','NUnsupportedOxygen', ...
     'NAboveUmin','NAtOrBelowUmin','NUnknownUmin','NBaseSolved', ...
     'NGAUKSSolved','NPairedBio','NFluxReactions', ...
     'NFluxMissingInBase','NFluxMissingInGAUKS'});

%% Per-reaction exchange-flux RMSE by model pair and condition group

pairKeys = unique(conditionResults(:, ...
    {'Species','BaseModel','GAUKSModel'}),'rows','stable');

groupNames = ["All","AboveUmin","AtOrBelowUmin"];
RReactionSummary = cell(0,14);

for p = 1:height(pairKeys)

    for gg = 1:numel(groupNames)

        groupName = groupNames(gg);

        maskPair = fluxResults.Species == pairKeys.Species(p) & ...
                   fluxResults.BaseModel == pairKeys.BaseModel(p) & ...
                   fluxResults.GAUKSModel == pairKeys.GAUKSModel(p);

        maskGroup = groupMask(fluxResults.UminGroup,groupName);
        supportedOxygen = fluxResults.OxAvail == "aerobic" | ...
                          fluxResults.OxAvail == "anaerobic";

        F = fluxResults(maskPair & maskGroup & supportedOxygen,:);

        if isempty(F)
            continue
        end

        reactions = unique(F.Reaction,'stable');

        for r = 1:numel(reactions)

            rxn = reactions(r);
            fr = F(F.Reaction == rxn,:);

            nMeasured = sum(isfinite(fr.FluxExperimental));
            paired = isfinite(fr.FluxExperimental) & ...
                     isfinite(fr.FluxBefore) & isfinite(fr.FluxAfter);
            nPaired = sum(paired);

            rmseB = NaN;
            rmseG = NaN;
            deltaRMSE = NaN;
            nrmseB = NaN;
            nrmseG = NaN;
            deltaNRMSE = NaN;
            scaleValues = fr.Scale(isfinite(fr.Scale) & fr.Scale > 0);
            scale = NaN;
            if ~isempty(scaleValues)
                scale = scaleValues(1);
            end
            if nPaired > 0
                y = fr.FluxExperimental(paired);
                predB = fr.FluxBefore(paired);
                predG = fr.FluxAfter(paired);
                rmseB = sqrt(mean((predB-y).^2));
                rmseG = sqrt(mean((predG-y).^2));
                deltaRMSE = rmseG-rmseB;
                if isfinite(scale) && scale > 0
                    nrmseB = rmseB/scale;
                    nrmseG = rmseG/scale;
                    deltaNRMSE = nrmseG-nrmseB;
                end
            end
            RReactionSummary(end+1,:) = { ...
                pairKeys.Species(p), pairKeys.BaseModel(p), ...
                pairKeys.GAUKSModel(p), groupName, rxn, ...
                nMeasured, nPaired, scale, rmseB, rmseG, ...
                deltaRMSE, nrmseB, nrmseG, deltaNRMSE};

        end
    end
end

reactionSummary = cell2table(RReactionSummary,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','Group','Reaction', ...
    'NMeasured','NPaired','Scale','RMSEBefore','RMSEAfter', ...
    'DeltaRMSE','NRMSEBefore','NRMSEAfter','DeltaNRMSE'});
reactionSummary = normalizeCellTable(reactionSummary, ...
    {'Species','BaseModel','GAUKSModel','Group','Reaction'}, ...
    {'NMeasured','NPaired','Scale','RMSEBefore','RMSEAfter', ...
     'DeltaRMSE','NRMSEBefore','NRMSEAfter','DeltaNRMSE'});

%% Per-model-pair summary: biomass, exchange macro-average, combined metric

RPairSummary = cell(0,25);
for p = 1:height(pairKeys)
    s = find(speciesNames == pairKeys.Species(p),1);
    bioScale = bioScales(s);
    for gg = 1:numel(groupNames)
        groupName = groupNames(gg);
        maskPairC = conditionResults.Species == pairKeys.Species(p) & ...
                    conditionResults.BaseModel == pairKeys.BaseModel(p) & ...
                    conditionResults.GAUKSModel == pairKeys.GAUKSModel(p);
        maskGroupC = groupMask(conditionResults.UminGroup,groupName);
        supportedC = conditionResults.OxAvail == "aerobic" | ...
                     conditionResults.OxAvail == "anaerobic";
        C = conditionResults(maskPairC & maskGroupC & supportedC,:);
        maskPairF = fluxResults.Species == pairKeys.Species(p) & ...
                    fluxResults.BaseModel == pairKeys.BaseModel(p) & ...
                    fluxResults.GAUKSModel == pairKeys.GAUKSModel(p);
        maskGroupF = groupMask(fluxResults.UminGroup,groupName);
        supportedF = fluxResults.OxAvail == "aerobic" | ...
                     fluxResults.OxAvail == "anaerobic";
        F = fluxResults(maskPairF & maskGroupF & supportedF,:);
        nConditions = height(C);
        pairedBio = isfinite(C.GrowthExperimental) & ...
                    isfinite(C.GrowthBefore) & isfinite(C.GrowthAfter);
        nPairedBio = sum(pairedBio);
        rmseBioB = NaN;
        rmseBioG = NaN;
        deltaBio = NaN;
        nrmseBioB = NaN;
        nrmseBioG = NaN;
        deltaBioNRMSE = NaN;
        if nPairedBio > 0
            y = C.GrowthExperimental(pairedBio);
            predB = C.GrowthBefore(pairedBio);
            predG = C.GrowthAfter(pairedBio);
            rmseBioB = sqrt(mean((predB-y).^2));
            rmseBioG = sqrt(mean((predG-y).^2));
            deltaBio = rmseBioG-rmseBioB;
            if isfinite(bioScale) && bioScale > 0
                nrmseBioB = rmseBioB/bioScale;
                nrmseBioG = rmseBioG/bioScale;
                deltaBioNRMSE = nrmseBioG-nrmseBioB;
            end
        end
        % Equal-weight macro-average over measured reactions.
        maskPairR = reactionSummary.Species == pairKeys.Species(p) & ...
                    reactionSummary.BaseModel == pairKeys.BaseModel(p) & ...
                    reactionSummary.GAUKSModel == pairKeys.GAUKSModel(p) & ...
                    reactionSummary.Group == groupName;
        RS = reactionSummary(maskPairR,:);
        useReaction = isfinite(RS.NRMSEBefore) & ...
                      isfinite(RS.NRMSEAfter) & RS.NPaired > 0;
        nReactionsUsed = sum(useReaction);
        exchangeB = NaN;
        exchangeG = NaN;
        deltaExchange = NaN;
        if nReactionsUsed > 0
            exchangeB = mean(RS.NRMSEBefore(useReaction));
            exchangeG = mean(RS.NRMSEAfter(useReaction));
            deltaExchange = exchangeG-exchangeB;
        end
        % Auxiliary combined score: biomass and exchange each contribute 50%.
        combinedB = NaN;
        combinedG = NaN;
        deltaCombined = NaN;
        if isfinite(nrmseBioB) && isfinite(nrmseBioG) && ...
                isfinite(exchangeB) && isfinite(exchangeG)
            combinedB = (nrmseBioB+exchangeB)/2;
            combinedG = (nrmseBioG+exchangeG)/2;
            deltaCombined = combinedG-combinedB;
        end
        ciBioLow = NaN; ciBioHigh = NaN;
        ciExchangeLow = NaN; ciExchangeHigh = NaN;
        ciCombinedLow = NaN; ciCombinedHigh = NaN;
        if nPairedBio >= 2
            [ciBioLow,ciBioHigh,ciExchangeLow,ciExchangeHigh, ...
             ciCombinedLow,ciCombinedHigh] = ...
                bootstrapPairMetrics(C,F,bioScale,nBootPair);
        end
        RPairSummary(end+1,:) = { ...
            pairKeys.Species(p), pairKeys.BaseModel(p), ...
            pairKeys.GAUKSModel(p), groupName, nConditions, ...
            nPairedBio, nReactionsUsed, ...
            rmseBioB, rmseBioG, deltaBio, ciBioLow, ciBioHigh, ...
            nrmseBioB, nrmseBioG, deltaBioNRMSE, ...
            exchangeB, exchangeG, deltaExchange, ...
            ciExchangeLow, ciExchangeHigh, ...
            combinedB, combinedG, deltaCombined, ...
            ciCombinedLow, ciCombinedHigh};
    end
end
pairSummary = cell2table(RPairSummary,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','Group','NConditions', ...
    'NPairedBio','NReactionsUsed','BioRMSEBefore','BioRMSEAfter', ...
    'DeltaBioRMSE','BioDeltaCI_Lower','BioDeltaCI_Upper', ...
    'BioNRMSEBefore','BioNRMSEAfter','DeltaBioNRMSE', ...
    'ExchangeNRMSEBefore','ExchangeNRMSEAfter','DeltaExchangeNRMSE', ...
    'ExchangeDeltaCI_Lower','ExchangeDeltaCI_Upper', ...
    'CombinedBefore','CombinedAfter','DeltaCombined', ...
    'CombinedDeltaCI_Lower','CombinedDeltaCI_Upper'});
pairSummary = normalizeCellTable(pairSummary, ...
    {'Species','BaseModel','GAUKSModel','Group'}, ...
    {'NConditions','NPairedBio','NReactionsUsed', ...
     'BioRMSEBefore','BioRMSEAfter','DeltaBioRMSE', ...
     'BioDeltaCI_Lower','BioDeltaCI_Upper', ...
     'BioNRMSEBefore','BioNRMSEAfter','DeltaBioNRMSE', ...
     'ExchangeNRMSEBefore','ExchangeNRMSEAfter','DeltaExchangeNRMSE', ...
     'ExchangeDeltaCI_Lower','ExchangeDeltaCI_Upper', ...
     'CombinedBefore','CombinedAfter','DeltaCombined', ...
     'CombinedDeltaCI_Lower','CombinedDeltaCI_Upper'});

%% Species-level summary across model pairs
metricNames = ["BioRMSE","ExchangeNRMSE","CombinedNRMSE"];
deltaColumns = {'DeltaBioRMSE','DeltaExchangeNRMSE','DeltaCombined'};
RSpeciesSummary = cell(0,12);
for s = 1:numel(speciesNames)
    for gg = 1:numel(groupNames)
        groupName = groupNames(gg);
        for mm = 1:numel(metricNames)
            deltaColumn = deltaColumns{mm};
            mask = pairSummary.Species == speciesNames(s) & ...
                   pairSummary.Group == groupName;
            valuesAll = pairSummary.(deltaColumn)(mask);
            values = valuesAll(isfinite(valuesAll));
            nPairsTotal = sum(mask);
            nPairsWithMetric = numel(values);
            meanDelta = NaN;
            medianDelta = NaN;
            q1Delta = NaN;
            q3Delta = NaN;
            improvedFraction = NaN;
            ciLow = NaN;
            ciHigh = NaN;
            if nPairsWithMetric > 0
                meanDelta = mean(values);
                medianDelta = median(values);
                q1Delta = simplePercentile(values,25);
                q3Delta = simplePercentile(values,75);
                improvedFraction = mean(values < 0);
                [ciLow,ciHigh] = bootstrapMeanCI(values,nBootSpecies);
            end
            RSpeciesSummary(end+1,:) = { ...
                speciesNames(s), groupName, metricNames(mm), ...
                nPairsTotal, nPairsWithMetric, meanDelta, medianDelta, ...
                q1Delta, q3Delta, improvedFraction, ciLow, ciHigh};

        end
    end
end
speciesSummary = cell2table(RSpeciesSummary,'VariableNames',{ ...
    'Species','Group','Metric','NModelPairsTotal', ...
    'NModelPairsWithResult','MeanDelta','MedianDelta', ...
    'Q1Delta','Q3Delta','ImprovedFraction', ...
    'BootstrapCI_Lower','BootstrapCI_Upper'});
speciesSummary = normalizeCellTable(speciesSummary, ...
    {'Species','Group','Metric'}, ...
    {'NModelPairsTotal','NModelPairsWithResult','MeanDelta', ...
     'MedianDelta','Q1Delta','Q3Delta','ImprovedFraction', ...
     'BootstrapCI_Lower','BootstrapCI_Upper'});

%% Export results
outDir = fullfile(findECOMAProot,'projects','Overall_Analysis');
if ~exist(outDir,'dir')
    mkdir(outDir);
end
outFile = fullfile(outDir,'GAUKS_Growth_Experiment_Full_Analysis.xlsx');
if exist(outFile,'file')
    delete(outFile);
end
writetable(conditionResults,outFile,'Sheet','ConditionPredictions');
writetable(fluxResults,outFile,'Sheet','FluxPredictions');
writetable(reactionSummary,outFile,'Sheet','ReactionSummary');
writetable(pairSummary,outFile,'Sheet','PairSummary');
writetable(speciesSummary,outFile,'Sheet','SpeciesSummary');
writetable(dataAudit,outFile,'Sheet','DataAudit');

disp("Growth experiment analysis complete.");
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
        if contains(ox(i),"anaerobic")
            ox(i) = "anaerobic";
        elseif contains(ox(i),"aerobic")
            ox(i) = "aerobic";
        end
    end
end
function scale = calculateScale(y)
    y = y(isfinite(y));
    if isempty(y)
        scale = NaN;
        return
    end
    scale = std(y,0);
    % If SD is zero/near-zero, use mean absolute magnitude as a fixed scale.
    if ~isfinite(scale) || scale <= eps(max(abs(y)))
        meanAbs = mean(abs(y));
        if isfinite(meanAbs) && meanAbs > 0
            scale = meanAbs;
        else
            scale = NaN;
        end
    end
end
function T = normalizeCellTable(T,stringVars,numericVars)
    for i = 1:numel(stringVars)
        v = stringVars{i};
        T.(v) = string(T.(v));
    end
    for i = 1:numel(numericVars)
        v = numericVars{i};
        x = T.(v);
        if iscell(x)
            T.(v) = cellfun(@(z) double(z),x);
        else
            T.(v) = double(x);
        end
    end
end
function mask = groupMask(groupColumn,groupName)
    if groupName == "All"
        mask = true(size(groupColumn));
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
    position = 1 + (numel(x)-1)*p/100;
    lo = floor(position);
    hi = ceil(position);
    if lo == hi
        value = x(lo);
    else
        value = x(lo) + (position-lo)*(x(hi)-x(lo));
    end
end
function [ciLow,ciHigh] = bootstrapMeanCI(values,nBoot)
    values = values(isfinite(values));
    if isempty(values)
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
    ciLow = simplePercentile(bootMeans,2.5);
    ciHigh = simplePercentile(bootMeans,97.5);
end
function [ciBioLow,ciBioHigh,ciExLow,ciExHigh,ciCombLow,ciCombHigh] = ...
    bootstrapPairMetrics(C,F,bioScale,nBoot)
    ciBioLow = NaN; ciBioHigh = NaN;
    ciExLow = NaN; ciExHigh = NaN;
    ciCombLow = NaN; ciCombHigh = NaN;
    % Resample only conditions with paired, observed biomass and both predictions.
    pairedBio = isfinite(C.GrowthExperimental) & ...
                isfinite(C.GrowthBefore) & isfinite(C.GrowthAfter);
    C = C(pairedBio,:);
    n = height(C);
    if n < 2
        return
    end
    yBio = C.GrowthExperimental;
    predBioB = C.GrowthBefore;
    predBioG = C.GrowthAfter;
    condIDs = C.ConditionID;
    reactions = unique(F.Reaction,'stable');
    nRxn = numel(reactions);
    % Align exchange-flux data into condition-by-reaction matrices.
    yFlux = NaN(n,nRxn);
    predFluxB = NaN(n,nRxn);
    predFluxG = NaN(n,nRxn);
    scales = NaN(1,nRxn);
    for r = 1:nRxn
        fr = F(F.Reaction == reactions(r),:);
        [tf,loc] = ismember(condIDs,fr.ConditionID);
        rows = find(tf);
        if ~isempty(rows)
            yFlux(rows,r) = fr.FluxExperimental(loc(rows));
            predFluxB(rows,r) = fr.FluxBefore(loc(rows));
            predFluxG(rows,r) = fr.FluxAfter(loc(rows));
        end
        sv = fr.Scale(isfinite(fr.Scale) & fr.Scale > 0);
        if ~isempty(sv)
            scales(r) = sv(1);
        end
    end
    bootBio = NaN(nBoot,1);
    bootExchange = NaN(nBoot,1);
    bootCombined = NaN(nBoot,1);
    for b = 1:nBoot
        idx = randi(n,n,1);
        % Biomass RMSE difference on paired conditions.
        rmseBioB = sqrt(mean((predBioB(idx)-yBio(idx)).^2));
        rmseBioG = sqrt(mean((predBioG(idx)-yBio(idx)).^2));
        bootBio(b) = rmseBioG-rmseBioB;
        nrmseBioB = NaN;
        nrmseBioG = NaN;
        if isfinite(bioScale) && bioScale > 0
            nrmseBioB = rmseBioB/bioScale;
            nrmseBioG = rmseBioG/bioScale;
        end
        exchangeB = [];
        exchangeG = [];
        % Each measured reaction contributes one equally weighted NRMSE.
        for r = 1:nRxn
            validFlux = isfinite(yFlux(idx,r)) & ...
                        isfinite(predFluxB(idx,r)) & ...
                        isfinite(predFluxG(idx,r));
            if ~any(validFlux) || ~isfinite(scales(r)) || scales(r) <= 0
                continue
            end
            yy = yFlux(idx(validFlux),r);
            ppB = predFluxB(idx(validFlux),r);
            ppG = predFluxG(idx(validFlux),r);
            rmseB = sqrt(mean((ppB-yy).^2));
            rmseG = sqrt(mean((ppG-yy).^2));
            exchangeB(end+1,1) = rmseB/scales(r);
            exchangeG(end+1,1) = rmseG/scales(r);
        end
        if ~isempty(exchangeB)
            exchangeNRMSEB = mean(exchangeB);
            exchangeNRMSEG = mean(exchangeG);
            bootExchange(b) = exchangeNRMSEG-exchangeNRMSEB;
            if isfinite(nrmseBioB) && isfinite(nrmseBioG)
                combinedB = (nrmseBioB+exchangeNRMSEB)/2;
                combinedG = (nrmseBioG+exchangeNRMSEG)/2;
                bootCombined(b) = combinedG-combinedB;
            end
        end
    end
    ciBioLow = simplePercentile(bootBio,2.5);
    ciBioHigh = simplePercentile(bootBio,97.5);
    ciExLow = simplePercentile(bootExchange,2.5);
    ciExHigh = simplePercentile(bootExchange,97.5);
    ciCombLow = simplePercentile(bootCombined,2.5);
    ciCombHigh = simplePercentile(bootCombined,97.5);
end
