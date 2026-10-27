%% Compare growth experiment predictions before and after GAUKS

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

%% Define species and models

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

ecGEM_bioRxns = {'biomass','r_2111','biomass'};

oxygenRxns = {'EX_o2_e','r_1992','EX_o2_e'};

%% Read experimental growth data

expSets = cell(1,3);

for s = 1:3
    expSets{s} = readtable( ...
        fullfile(modelLocations{s}, ...
        'data','BayesianGrowthRates.tsv'), ...
        'FileType','text', ...
        'Delimiter','\t', ...
        'VariableNamingRule','preserve');
end

%% Parameters

uptakeTol = 1e-6;
nBoot = 2000;

% Condition-level results
R = cell(0,20);

%% Evaluate each species

for s = 1:numel(speciesNames)

    M = modelSets{s};
    T = expSets{s};

    % Experimental measurements
    target = double(T.(ecGEM_bioRxns{s}));
    uptakeAll = double(T.Uptake);

    conditionID = (1:height(T))';

    valid = isfinite(target) & ...
            isfinite(uptakeAll) & ...
            uptakeAll < 0;

    T = T(valid,:);
    target = target(valid);
    uptakeAll = uptakeAll(valid);
    conditionID = conditionID(valid);

    % Carbon substrates represented in the experimental data
    allCarbon = unique(string(T.Substrate),'stable');

    allCarbon = allCarbon(~ismissing(allCarbon) & ...
                          allCarbon ~= "" & ...
                          allCarbon ~= "NaN");

    % Find GAUKS models and their corresponding base models
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

        % Set biomass objectives
        ibB = findRxnIDs(mB,ecGEM_bioRxns{s});
        ibG = findRxnIDs(mG,ecGEM_bioRxns{s});

        if ibB == 0 || ibG == 0
            continue
        end

        mB.c(:) = 0;
        mG.c(:) = 0;

        mB.c(ibB) = 1;
        mG.c(ibG) = 1;

        % Identify substrate exchange reactions
        idsB = findRxnIDs(mB,cellstr(allCarbon));
        idsG = findRxnIDs(mG,cellstr(allCarbon));

        % Extract GAUKS parameters
        params = extractSluiceParams(mG,cellstr(allCarbon));

        %% Evaluate each experimental condition

        for i = 1:height(T)

            ex = string(T.Substrate(i));
            ox = lower(strtrim(string(T.OxAvail(i))));
            media = string(T.Media(i));

            muExp = target(i);
            uptakeExp = uptakeAll(i);
            uptakeMagnitude = abs(uptakeExp);

            % Identify substrate-specific Umin
            k = find(allCarbon == ex,1);

            umin = NaN;

            if ~isempty(k)
                umin = params.umin(k);
            end

            % Classify experimental uptake relative to Umin
            if ~isfinite(umin)
                uminGroup = "UnknownUmin";
            elseif uptakeMagnitude > umin + uptakeTol
                uminGroup = "AboveUmin";
            elseif uptakeMagnitude < umin - uptakeTol
                uminGroup = "BelowUmin";
            else
                uminGroup = "AtUmin";
            end

            % Initialize predictions
            muB = NaN;
            muG = NaN;

            statusB = "NotSolved";
            statusG = "NotSolved";

            % Identify substrate reactions
            jB = findRxnIDs(mB,char(ex));
            jG = findRxnIDs(mG,char(ex));

            iBasal = findRxnIDs(mG,char(ex+"_basal"));
            iExtended = findRxnIDs(mG,char(ex+"_extended"));

            % This script evaluates aerobic and anaerobic conditions.
            % "limited" oxygen conditions are excluded because this
            % dataset does not provide a uniform oxygen uptake bound.
            if ox ~= "aerobic" && ox ~= "anaerobic"

                statusB = "Excluded_OxygenCondition";
                statusG = "Excluded_OxygenCondition";

            else

                %% Base model

                if jB == 0

                    statusB = "MissingSubstrateRxn";

                else

                    b = mB;

                    % Close uptake of other listed carbon substrates
                    idxB = idsB(idsB > 0);

                    b.lb(idxB) = 0;
                    b.ub(idxB) = max(0,b.ub(idxB));

                    % Fix substrate uptake to the experimental value
                    b.lb(jB) = uptakeExp;
                    b.ub(jB) = uptakeExp;

                    % Apply anaerobic conditions
                    if ox == "anaerobic"
                        b = anaerobicModel( ...
                            b,orgNames{s},oxygenRxns{s});
                    end

                    % Maximize growth
                    sb = solveLP(b,'max');

                    if ~isempty(sb.x) && sb.stat == 1
                        muB = sb.x(ibB);
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

                    % Close uptake of other listed carbon substrates
                    idxG = idsG(idsG > 0);

                    g.lb(idxG) = 0;
                    g.ub(idxG) = max(0,g.ub(idxG));

                    % Fix the same substrate uptake as the base model
                    g.lb(jG) = uptakeExp;
                    g.ub(jG) = uptakeExp;

                    % Apply anaerobic conditions
                    if ox == "anaerobic"
                        g = anaerobicModel( ...
                            g,orgNames{s},oxygenRxns{s});
                    end

                    % Maximize growth
                    sg = solveLP(g,'max');

                    if ~isempty(sg.x) && sg.stat == 1
                        muG = sg.x(ibG);
                        statusG = "Solved";
                    else
                        statusG = "SolveFailed";
                    end

                end
            end

            %% Calculate errors

            errorB = muB - muExp;
            errorG = muG - muExp;

            absErrorB = abs(errorB);
            absErrorG = abs(errorG);

            % Only aerobic and anaerobic conditions are eligible
            eligible = ox == "aerobic" || ox == "anaerobic";

            %% Save condition-level results

            R(end+1,:) = { ...
                speciesNames(s), ...
                string(baseName), ...
                string(gauksName), ...
                conditionID(i), ...
                ex, ...
                ox, ...
                media, ...
                uptakeExp, ...
                uptakeMagnitude, ...
                muExp, ...
                umin, ...
                uminGroup, ...
                muB, ...
                muG, ...
                errorB, ...
                errorG, ...
                absErrorB, ...
                absErrorG, ...
                statusB, ...
                statusG};

        end
    end
end

%% Convert condition-level results to table

growthResults = cell2table(R,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','ConditionID', ...
    'Substrate','OxAvail','Media', ...
    'UptakeExperimental','UptakeMagnitude', ...
    'GrowthExperimental','Umin','UminGroup', ...
    'GrowthBefore','GrowthAfter', ...
    'ErrorBefore','ErrorAfter', ...
    'AbsErrorBefore','AbsErrorAfter', ...
    'BaseStatus','GAUKSStatus'});

%% Summarize RMSE by species, model pair and Umin group

pairKeys = unique(growthResults(:, ...
    {'Species','BaseModel','GAUKSModel'}), ...
    'rows','stable');

groupNames = ["All","BelowUmin","AboveUmin"];

S = cell(0,15);

for p = 1:height(pairKeys)

    maskPair = ...
        growthResults.Species == pairKeys.Species(p) & ...
        growthResults.BaseModel == pairKeys.BaseModel(p) & ...
        growthResults.GAUKSModel == pairKeys.GAUKSModel(p);

    for gg = 1:numel(groupNames)

        groupName = groupNames(gg);

        if groupName == "All"
            maskGroup = true(height(growthResults),1);
        else
            maskGroup = growthResults.UminGroup == groupName;
        end

        mask = maskPair & maskGroup & ...
            isfinite(growthResults.GrowthExperimental);

        % Exclude oxygen conditions without a defined constraint
        % from the primary comparison
        ox = growthResults.OxAvail;
        mask = mask & (ox == "aerobic" | ox == "anaerobic");

        nConditions = sum(mask);

        maskB = mask & isfinite(growthResults.GrowthBefore);
        maskG = mask & isfinite(growthResults.GrowthAfter);

        maskPaired = mask & ...
            isfinite(growthResults.GrowthBefore) & ...
            isfinite(growthResults.GrowthAfter);

        nBefore = sum(maskB);
        nAfter = sum(maskG);
        nPaired = sum(maskPaired);

        rmseB = NaN;
        rmseG = NaN;
        deltaRMSE = NaN;
        ciLow = NaN;
        ciHigh = NaN;

        if nPaired > 0

            y = growthResults.GrowthExperimental(maskPaired);
            predB = growthResults.GrowthBefore(maskPaired);
            predG = growthResults.GrowthAfter(maskPaired);

            % Both RMSE values use the same paired conditions
            rmseB = sqrt(mean((predB-y).^2));
            rmseG = sqrt(mean((predG-y).^2));

            % Negative values indicate lower RMSE after GAUKS
            deltaRMSE = rmseG-rmseB;

            if nPaired >= 2
                [ciLow,ciHigh] = pairedBootstrapDeltaRMSE( ...
                    y,predB,predG,nBoot);
            end

        end

        coverageB = NaN;
        coverageG = NaN;

        if nConditions > 0
            coverageB = nBefore/nConditions;
            coverageG = nAfter/nConditions;
        end

        S(end+1,:) = { ...
            pairKeys.Species(p), ...
            pairKeys.BaseModel(p), ...
            pairKeys.GAUKSModel(p), ...
            groupName, ...
            nConditions, ...
            nBefore, ...
            nAfter, ...
            nPaired, ...
            rmseB, ...
            rmseG, ...
            deltaRMSE, ...
            ciLow, ...
            ciHigh, ...
            coverageB, ...
            coverageG};

    end
end

%% Convert summary results to table

growthSummary = cell2table(S,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','Group', ...
    'NConditions','NBeforeSolved','NAfterSolved','NPaired', ...
    'RMSEBefore','RMSEAfter','DeltaRMSE', ...
    'BootstrapCI_Lower','BootstrapCI_Upper', ...
    'CoverageBefore','CoverageAfter'});

%% Export results

outFile = fullfile(findECOMAProot, ...
    'projects','Overall_Analysis', ...
    'GAUKS_Growth_Experiment_Comparison.xlsx');

writetable(growthResults,outFile,'Sheet','GrowthPredictions');
writetable(growthSummary,outFile,'Sheet','GrowthRMSE');

%% Paired bootstrap function

function [ciLow,ciHigh] = pairedBootstrapDeltaRMSE( ...
    y,predB,predG,nBoot)

    n = numel(y);
    bootDelta = zeros(nBoot,1);

    for b = 1:nBoot

        % Resample identical condition indices for both models
        idx = randi(n,n,1);

        rmseB = sqrt(mean((predB(idx)-y(idx)).^2));
        rmseG = sqrt(mean((predG(idx)-y(idx)).^2));

        bootDelta(b) = rmseG-rmseB;

    end

    ci = prctile(bootDelta,[2.5 97.5]);

    ciLow = ci(1);
    ciHigh = ci(2);

end