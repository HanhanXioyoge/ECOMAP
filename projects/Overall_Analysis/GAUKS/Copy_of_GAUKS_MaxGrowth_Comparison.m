% Compare_before_and_after_adding_GAUKS
%% load models
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

%% Read data
speciesNames = ["eciML1515","ecYeast","eciCW773"];
orgNames = {'Escherichia coli', ...
            'Saccharomyces cerevisiae', ...
            'Corynebacterium glutamicum'};
modelSets = {eciML1515_models,ecYeast_models,eciCW773_models};
modelLocations = {eciML1515_kcatRepo_Location,...
                  ecYeast_kcatRepo_Location,...
                  eciCW773_kcatRepo_Location};
ecGEM_bioRxns = {'biomass','r_2111','biomass'};
cSources = {'EX_glc__D_e', 'r_1714', 'EX_glc_e'};
oxygenRxns = {'EX_o2_e','r_1992','EX_o2_e'};
tableSets = cell(1,3);
expSets = cell(1,3);
for s = 1:3
    tableSets{s} = readtable( ...
        fullfile(modelLocations{s},'data','UnconstrainedMaxGrowth.tsv'), ...
        'FileType','text','Delimiter','\t','VariableNamingRule','preserve');
    expSets{s} = readtable( ...
        fullfile(modelLocations{s},'data','BayesianGrowthRates.tsv'), ...
        'FileType','text','Delimiter','\t','VariableNamingRule','preserve');
end
uptakeCap = 1000;
growthTol = 1e-6;
paramTol = 1e-9;
R = cell(0,17);
for s = 1:numel(speciesNames)
    M = modelSets{s};
    T = tableSets{s};
    target = double(T.(ecGEM_bioRxns{s}));
    valid = isfinite(target) & target > 0;
    T = T(valid,:);
    target = target(valid);
    allCarbon = unique(string(T.Substrate),'stable');
    allCarbon = allCarbon(~ismissing(allCarbon) & ...
        allCarbon ~= "NaN" & allCarbon ~= "");
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
        ibB = findRxnIDs(mB,ecGEM_bioRxns{s});
        ibG = findRxnIDs(mG,ecGEM_bioRxns{s});
        if ibB == 0 || ibG == 0
            continue
        end
        % Set biomass objective
        mB.c(:) = 0;
        mG.c(:) = 0;
        mB.c(ibB) = 1;
        mG.c(ibG) = 1;
        % Identify carbon exchange reactions
        idsB = findRxnIDs(mB,cellstr(allCarbon));
        idsG = findRxnIDs(mG,cellstr(allCarbon));
        % Extract GAUKS parameters
        params = extractSluiceParams(mG,cellstr(allCarbon));
        %% Maximum growth and relative error
        for i = 1:height(T)
            ex = string(T.Substrate(i));
            ox = string(T.OxAvail(i));
            targetMu = target(i);
            k = find(allCarbon == ex,1);
            if isempty(k)
                continue
            end
            jB = idsB(k);
            jG = idsG(k);
            iBasal = findRxnIDs(mG,char(ex+"_basal"));
            iExtended = findRxnIDs(mG,char(ex+"_extended"));
            if jB == 0 || jG == 0 || ...
                    iBasal == 0 || iExtended == 0
                continue
            end
            umin = params.umin(k);
            xi = params.xi(k);
            % Copy models
            b = mB;
            g = mG;
            % Close uptake of all other carbon sources
            b.lb(idsB(idsB>0)) = 0;
            b.ub(idsB(idsB>0)) = max(0,b.ub(idsB(idsB>0)));
            g.lb(idsG(idsG>0)) = 0;
            g.ub(idsG(idsG>0)) = max(0,g.ub(idsG(idsG>0)));
            % Allow uptake of the tested substrate only
            b.lb(jB) = -uptakeCap;
            b.ub(jB) = 0;
            g.lb(jG) = -uptakeCap;
            g.ub(jG) = 0;
            % Apply anaerobic conditions when required
            if ox == "anaerobic"
                b = anaerobicModel(b,orgNames{s},oxygenRxns{s});
                g = anaerobicModel(g,orgNames{s},oxygenRxns{s});
            end
            % Maximize growth
            sb = solveLP(b,'max');
            sg = solveLP(g,'max');
            muB = NaN;
            muG = NaN;
            qB = NaN;
            qG = NaN;
            if ~isempty(sb.x) && sb.stat == 1
                muB = sb.x(ibB);
                qB = max(0,-sb.x(jB));
            end
            if ~isempty(sg.x) && sg.stat == 1
                muG = sg.x(ibG);
                qG = max(0,-sg.x(jG));
            end
            % Determine whether GAUKS was activated
            searchTriggered = isfinite(muB) && ...
                isfinite(targetMu) && ...
                muB > targetMu + growthTol;
            if ~isfinite(umin) || umin < 0
                status = "InvalidUmin";
            elseif umin >= uptakeCap-paramTol
                status = "Inactive_Umin1000";
            elseif ~isfinite(xi)
                status = "InvalidXi";
            elseif xi <= paramTol
                status = "Inactive_XiZero";
            elseif ~isfinite(muB)
                status = "BaseSolveFailed";
            elseif ~isfinite(targetMu) || targetMu <= 0
                status = "InvalidGrowthReference";
            elseif ~isfinite(muG)
                status = "GAUKSSolveFailed";
            elseif ~searchTriggered
                status = "Inactive_NoOverprediction";
            else
                status = "Active";
            end
            % Growth rate change
            growthChange = muG - muB;
            growthChangePct = NaN;
            if isfinite(muB) && abs(muB) > growthTol
                growthChangePct = 100*growthChange/muB;
            end
            % Relative growth errors: omega
            omegaB = NaN;
            omegaG = NaN;
            deltaOmega = NaN;
            if isfinite(targetMu) && targetMu > 0
                if isfinite(muB)
                    omegaB = abs(muB-targetMu)/targetMu;
                end

                if isfinite(muG)
                    omegaG = abs(muG-targetMu)/targetMu;
                end

                if isfinite(omegaB) && isfinite(omegaG)
                    deltaOmega = omegaG-omegaB;
                end
            end
            % Store results
            R(end+1,:) = { ...
                speciesNames(s), ...
                string(baseName), ...
                string(gauksName), ...
                ex, ...
                ox, ...
                targetMu, ...
                muB, ...
                muG, ...
                growthChange, ...
                growthChangePct, ...
                qB, ...
                qG, ...
                qG-qB, ...
                omegaB, ...
                omegaG, ...
                deltaOmega, ...
                status};

        end
    end
end
%% Convert results to table
growthResults = cell2table(R,'VariableNames',{ ...
    'Species','BaseModel','GAUKSModel','Substrate','OxAvail', ...
    'GrowthReference','GrowthBefore','GrowthAfter', ...
    'GrowthChange','GrowthChangePercent', ...
    'UptakeBefore','UptakeAfter','UptakeChange', ...
    'OmegaBefore','OmegaAfter','DeltaOmega','GAUKSStatus'});
%% Export results
outFile = fullfile(findECOMAProot,'projects', 'Overall_Analysis',...
    'GAUKS_MaxGrowth_Comparison.xlsx');
writetable(growthResults,outFile,'Sheet','MaxGrowth');
