%% Physiological uptake diagnosis of the unmodified base ecGEMs
% Run setup from the repository root first. No GAUKS or calibration is run.
%% Load models (same three-species configuration as the comparison scripts)
eciML1515_Location = fullfile(findECOMAProot,'projects','eciML1515');
ecYeast_Location = fullfile(findECOMAProot,'projects','ecYeast');
eciCW773_Location = fullfile(findECOMAProot,'projects','eciCW773');
load(fullfile(eciML1515_Location,'models','eciML1515_models.mat'))
load(fullfile(ecYeast_Location,'models','ecYeast_models.mat'))
load(fullfile(eciCW773_Location,'models','eciCW773_models.mat'))

speciesNames = ["eciML1515","ecYeast","eciCW773"];
modelSets = {eciML1515_models,ecYeast_models,eciCW773_models};
modelLocations = {fullfile(eciML1515_Location,'calibration'), ...
    fullfile(ecYeast_Location,'calibration'), ...
    fullfile(eciCW773_Location,'calibration')};
ecGEM_bioRxns = {'biomass','r_2111','biomass'};
oxygenRxns = {'EX_o2_e','r_1992','EX_o2_e'};
ammoniumRxns = {{'EX_nh4_e','EX_nh3_e'}, ...
    {'r_1654'}, {'EX_nh4_e','EX_nh3_e'}};

%% Read the original condition tables; keep their row IDs and row order
nSpecies = numel(speciesNames);
expSets = cell(1,nSpecies);
tableSets = cell(1,nSpecies);
for s = 1:nSpecies
    tableSets{s} = readtable(fullfile(modelLocations{s},'data', ...
        'UnconstrainedMaxGrowth.tsv'), 'FileType','text', ...
        'Delimiter','\t','VariableNamingRule','preserve');
    expSets{s} = readtable(fullfile(modelLocations{s},'data', ...
        'BayesianGrowthRates.tsv'), 'FileType','text', ...
        'Delimiter','\t','VariableNamingRule','preserve');
end

pairedNames = {'Species','Model','Substrate','OxAvail','ConditionID', ...
    'GrowthExp','UptakeExp','UptakeMin','GrowthFixedUptake', ...
    'GrowthError','OxygenUptakePred_A','OxygenUptakePred_B', ...
    'NH3UptakePred_A','NH3UptakePred_B','ByproductFluxes', ...
    'LPStatusA','LPStatusB'};
maxNames = {'Species','Model','Substrate','OxAvail','GrowthMax', ...
    'UptakeAtGrowthMax','OxygenUptakeAtGrowthMax', ...
    'NH3UptakeAtGrowthMax','LPStatusC'};
pairedRows = cell(0,numel(pairedNames));
maxRows = cell(0,numel(maxNames));

for s = 1:nSpecies
    M = modelSets{s};
    E = expSets{s};
    U = tableSets{s};
    bioRxn = ecGEM_bioRxns{s};
    % GEM is not enzyme constrained. Every remaining field without _GAUKS
    % is a base ecGEM, including Init and Uncalibrated where present.
    modelNames = fieldnames(M);
    modelNames = modelNames(~endsWith(modelNames,'_GAUKS') & ...
        ~strcmp(modelNames,'GEM'));
    % Existing project convention: each MIN row supplies one of these
    % alternative carbon sources. Leave all other medium bounds intact.
    carbonRxns = unique([string(E.Substrate);string(U.Substrate)],'stable');
    carbonRxns = carbonRxns(~ismissing(carbonRxns) & carbonRxns ~= "" & ...
        lower(carbonRxns) ~= "nan");

    for p = 1:numel(modelNames)
        modelName = modelNames{p};
        base = M.(modelName);
        ib = findRxnIDs(base,bioRxn);
        io = findRxnIDs(base,oxygenRxns{s});
        in = findRxnIDs(base,ammoniumRxns{s});
        candidateIds = findRxnIDs(base,cellstr(carbonRxns));

        for i = 1:height(E)
            target = string(E.Substrate(i));
            ox = lower(strtrim(string(E.OxAvail(i))));
            media = string(E.Media(i));
            muExp = numericValue(E{i,bioRxn});
            uptakeRaw = numericValue(E{i,'Uptake'});
            uptakeExp = abs(uptakeRaw);
            if ~isfinite(uptakeRaw), uptakeExp = NaN; end
            conditionID = string(E{i,1});
            jt = findRxnIDs(base,char(target));
            byproducts = measuredProducts(E,i,base,jt,ib,io,in);
            byproductText = productDescription(byproducts);
            uptakeMin = NaN; growthB = NaN; errorB = NaN;
            oA = NaN; oB = NaN; nA = NaN; nB = NaN;

            commonStatus = validateCondition(base,ib,jt,io,ox, ...
                isfinite(muExp) && muExp >= 0 && ...
                isfinite(uptakeExp) && uptakeExp > 0);
            if commonStatus ~= "OK"
                statusA = commonStatus; statusB = commonStatus;
            else
                % Test A: independent copy, fixed growth and measured
                % secretion only, minimum nonnegative target uptake.
                a = conditionModel(base,io,ox,media,candidateIds,jt);
                if muExp < a.lb(ib)-1e-9 || muExp > a.ub(ib)+1e-9
                    statusA = "ExperimentalGrowthOutsideBounds";
                else
                    a.lb(ib) = muExp; a.ub(ib) = muExp;
                    [a,productStatus] = fixMeasuredProducts(a,byproducts);
                    if productStatus ~= "OK"
                        statusA = productStatus;
                    else
                        direction = exchangeUptakeSign(a,jt);
                        if direction < 0
                            a.ub(jt) = min(a.ub(jt),0);
                        else
                            a.lb(jt) = max(a.lb(jt),0);
                        end
                        if a.lb(jt) > a.ub(jt)
                            statusA = "TargetUptakeBoundsConflict";
                        else
                            a.c(:) = 0; a.c(jt) = direction;
                            [vA,statusA] = solveBaseLP(a,'min');
                            if statusA == "Optimal"
                                uptakeMin = positiveUptake(a,jt,vA);
                                oA = totalUptake(a,io,vA);
                                nA = totalUptake(a,in,vA);
                            end
                        end
                    end
                end

                % Test B: independent copy; all products and biomass free.
                b = conditionModel(base,io,ox,media,candidateIds,jt);
                direction = exchangeUptakeSign(b,jt);
                fixedFlux = direction * uptakeExp;
                if fixedFlux < b.lb(jt)-1e-9 || fixedFlux > b.ub(jt)+1e-9
                    statusB = "ExperimentalUptakeOutsideBounds";
                else
                    b.lb(jt) = fixedFlux; b.ub(jt) = fixedFlux;
                    b.c(:) = 0; b.c(ib) = 1;
                    [vB,statusB] = solveBaseLP(b,'max');
                    if statusB == "Optimal"
                        growthB = vB(ib);
                        errorB = growthB-muExp;
                        oB = totalUptake(b,io,vB);
                        nB = totalUptake(b,in,vB);
                    end
                end
            end
            pairedRows(end+1,:) = {speciesNames(s),string(modelName), ...
                target,ox,conditionID,muExp,uptakeExp,uptakeMin,growthB, ...
                errorB,oA,oB,nA,nB,byproductText,statusA,statusB}; %#ok<SAGROW>
        end

        for i = 1:height(U)
            target = string(U.Substrate(i));
            ox = lower(strtrim(string(U.OxAvail(i))));
            media = string(U.Media(i));
            jt = findRxnIDs(base,char(target));
            growthMax = NaN; uptakeC = NaN; oC = NaN; nC = NaN;
            statusC = validateCondition(base,ib,jt,io,ox,true);
            if statusC == "OK"
                % Test C: target keeps its original exchange bounds.
                c = conditionModel(base,io,ox,media,candidateIds,jt);
                c.c(:) = 0; c.c(ib) = 1;
                [vC,statusC] = solveBaseLP(c,'max');
                if statusC == "Optimal"
                    growthMax = vC(ib);
                    uptakeC = positiveUptake(c,jt,vC);
                    oC = totalUptake(c,io,vC);
                    nC = totalUptake(c,in,vC);
                end
            end
            maxRows(end+1,:) = {speciesNames(s),string(modelName), ...
                target,ox,growthMax,uptakeC,oC,nC,statusC}; %#ok<SAGROW>
        end
    end
end

pairedResults = cell2table(pairedRows,'VariableNames',pairedNames);
maxGrowthResults = cell2table(maxRows,'VariableNames',maxNames);
outFile = fullfile(findECOMAProot,'projects','Overall_Analysis', ...
    'BaseModel_PhysiologicalUptake.xlsx');
writetable(pairedResults,outFile,'Sheet','PairedGrowthUptake');
writetable(maxGrowthResults,outFile,'Sheet','UnconstrainedMaxGrowth');
fprintf('Saved %d paired rows and %d maximum-growth rows to %s\n', ...
    height(pairedResults),height(maxGrowthResults),outFile);
disp(groupsummary(pairedResults,{'Species','LPStatusA','LPStatusB'}));
disp(groupsummary(maxGrowthResults,{'Species','LPStatusC'}));

%% Local helpers
function value = numericValue(x)
    if iscell(x), x = x{1}; end
    if isnumeric(x) || islogical(x)
        value = double(x);
    else
        value = str2double(string(x));
    end
    if isempty(value), value = NaN; end
end

function status = validateCondition(model,ib,jt,io,ox,validInput)
    status = "OK";
    if ~validInput || ismissing(ox) || ox == ""
        status = "MissingExperimentalInput";
    elseif ib == 0
        status = "MissingBiomassReaction";
    elseif jt == 0
        status = "MissingSubstrateReaction";
    elseif io == 0
        status = "MissingOxygenReaction";
    elseif ox == "limited"
        % No quantitative oxygen cap is supplied for this condition.
        status = "UnknownLimitedOxygenBound";
    elseif ox ~= "aerobic" && ox ~= "anaerobic"
        status = "UnknownOxygenCondition";
    elseif ~isfinite(exchangeUptakeSign(model,jt)) || ...
            ~isfinite(exchangeUptakeSign(model,io))
        status = "UnknownExchangeDirection";
    end
end

function model = conditionModel(model,io,ox,media,carbonIds,jt)
    % MIN data use one listed carbon source per condition. Close uptake of
    % alternative listed sources, retaining their secretion direction and
    % every other original medium, enzyme, and protein-pool constraint.
    if strcmpi(strtrim(media),'MIN')
        alternatives = setdiff(carbonIds(carbonIds > 0),jt);
        for k = 1:numel(alternatives)
            j = alternatives(k);
            d = exchangeUptakeSign(model,j);
            if d < 0, model.lb(j) = max(model.lb(j),0); end
            if d > 0, model.ub(j) = min(model.ub(j),0); end
        end
    end
    if ox == "anaerobic"
        d = exchangeUptakeSign(model,io);
        if d < 0, model.lb(io) = max(model.lb(io),0); end
        if d > 0, model.ub(io) = min(model.ub(io),0); end
    end
end

function d = exchangeUptakeSign(model,j)
    d = NaN;
    if numel(j) ~= 1 || j == 0, return; end
    coefficient = nonzeros(model.S(:,j));
    if isscalar(coefficient) && coefficient ~= 0
        d = sign(coefficient);
    end
end

function products = measuredProducts(T,i,model,jt,ib,io,in)
    products = cell(0,3); % reaction ID, model index, secretion magnitude
    skip = [jt,ib,io,in(:)'];
    vars = T.Properties.VariableNames;
    for k = 1:numel(vars)
        name = vars{k};
        if any(strcmp(name,{'Substrate','Uptake','OxAvail','Media'}))
            continue
        end
        j = findRxnIDs(model,name);
        if j == 0 || any(j == skip), continue; end
        raw = numericValue(T{i,name});
        if ~isfinite(raw), continue; end
        % Negative observations are uptake, not product secretion.
        d = exchangeUptakeSign(model,j);
        if ~isfinite(d) || raw < 0, continue; end
        products(end+1,:) = {name,j,raw}; %#ok<AGROW>
    end
end

function description = productDescription(products)
    description = "";
    for k = 1:size(products,1)
        if k > 1, description = description + "; "; end
        description = description + string(products{k,1}) + "=" + ...
            string(products{k,3});
    end
end

function [model,status] = fixMeasuredProducts(model,products)
    status = "OK";
    for k = 1:size(products,1)
        j = products{k,2};
        flux = -exchangeUptakeSign(model,j) * products{k,3};
        if flux < model.lb(j)-1e-9 || flux > model.ub(j)+1e-9
            status = "ProductFluxOutsideBounds:" + string(products{k,1});
            return
        end
        model.lb(j) = flux; model.ub(j) = flux;
    end
end

function [v,status] = solveBaseLP(model,sense)
    v = [];
    try
        % RAVEN solveLP always maximizes model.c. Its second argument is
        % minFlux (secondary flux minimization), not the objective sense.
        if strcmp(sense,'min'), model.c = -model.c; end
        sol = solveLP(model,0);
        if ~isstruct(sol) || ~isfield(sol,'stat') || isempty(sol.stat)
            status = "SolverReturnedNoStatus";
        elseif sol.stat == 1 && isfield(sol,'x') && ...
                numel(sol.x) == numel(model.rxns) && ...
                all(isfinite(sol.x))
            v = sol.x;
            status = "Optimal";
        elseif sol.stat == 0
            status = "FeasibleNotOptimal";
        elseif sol.stat == -1
            % solveLP collapses infeasible and unbounded to stat=-1.
            % A zero-objective feasibility solve distinguishes them.
            feasibilityModel = model;
            feasibilityModel.c(:) = 0;
            feasible = solveLP(feasibilityModel,0);
            if feasible.stat >= 0 && ~isempty(feasible.x)
                status = "Unbounded";
            else
                status = "Infeasible";
            end
        else
            status = "SolverStat_" + string(sol.stat);
        end
    catch ME
        status = "SolverError:" + string(ME.identifier);
    end
end

function q = positiveUptake(model,j,v)
    q = NaN;
    if j > 0 && isfinite(exchangeUptakeSign(model,j))
        q = max(0,exchangeUptakeSign(model,j)*v(j));
    end
end

function q = totalUptake(model,ids,v)
    q = NaN;
    ids = ids(ids > 0);
    if isempty(ids), return; end
    q = 0;
    for k = 1:numel(ids)
        value = positiveUptake(model,ids(k),v);
        if ~isfinite(value), q = NaN; return; end
        q = q + value;
    end
end
