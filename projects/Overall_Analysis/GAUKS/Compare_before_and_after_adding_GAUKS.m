% Compare before and after adding GAUKS
% Load models
eciML1515_Location = fullfile(findECOMAProot, 'projects', 'eciML1515'); 
eciML1515_models_Location = fullfile(eciML1515_Location, 'models'); 
eciML1515_kcatRepo_Location = fullfile(eciML1515_Location, 'calibration');
load(fullfile(eciML1515_models_Location, 'eciML1515_models_part1.mat'))
load(fullfile(eciML1515_models_Location, 'eciML1515_models_part2.mat'))
eciML1515_models = cell2struct([struct2cell(models_part1); struct2cell(models_part2)], ...
    [fieldnames(models_part1); fieldnames(models_part2)], 1);

ecYeast_Location = fullfile(findECOMAProot, 'projects', 'ecYeast'); 
ecYeast_models_Location = fullfile(ecYeast_Location, 'models'); 
ecYeast_kcatRepo_Location = fullfile(ecYeast_Location, 'calibration');
load(fullfile(ecYeast_models_Location, 'ecYeast_models.mat'))

eciCW773_Location = fullfile(findECOMAProot, 'projects', 'eciCW773'); 
eciCW773_models_Location = fullfile(eciCW773_Location, 'models'); 
eciCW773_kcatRepo_Location = fullfile(eciCW773_Location, 'calibration');
load(fullfile(eciCW773_models_Location, 'eciCW773_models.mat'))
%  Compare before and after adding GAUKS
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
% 13C data: currently available for eciML1515 only
C13Sets = cell(1,3);
C13Sets{1} = readtable(fullfile(eciML1515_kcatRepo_Location,'data','13CFluxdata.tsv'),'FileType','text','Delimiter','\t','VariableNamingRule','preserve');
C13Sets{2} = readtable(fullfile(ecYeast_kcatRepo_Location,'data','13CFluxdata.tsv'),'FileType','text','Delimiter','\t','VariableNamingRule','preserve');

%% Parameters
uptakeCap = 1000;
paramTol = 1e-8;
growthTol = 1e-7;
ratios = [0.25 0.50 0.75 0.90 0.99 1 1.01 1.10 1.25 1.50 2];
growthFractions = [0.25 0.50 0.75 0.90];
% All results are accumulated here and written to one Excel workbook.
R = cell(0,17);  % Maximum growth
P = cell(0,15);  % Parameter audit
A = cell(0,18);  % Threshold scan
B = cell(0,21);  % Experimental growth
C = cell(0,21);  % Uptake FVA
D = cell(0,13);  % Physiological feasibility
E = cell(0,16);  % 13C internal flux rmse
Z = cell(0,9);   % Zero-penalty control

%% Run all GAUKS analyses
for s = 1:3
    M = modelSets{s};
    T = tableSets{s};
    X = expSets{s};
    if s <= 2
        C13 = C13Sets{s};
        c13CondCols = find(startsWith(string(C13.Properties.VariableNames),"Cond"));
        c13ConditionNames = C13.Properties.VariableNames(c13CondCols);
        c13ConstraintRows = find(string(C13.Type)=="constraint");
        c13FluxRows = find(string(C13.Type)=="flux");
    else
        C13 = table();
        c13CondCols = [];
        c13ConditionNames = {};
        c13ConstraintRows = [];
        c13FluxRows = [];
    end
    T.Substrate = strtrim(erase(string(T.Substrate),"'"));
    X.Substrate = strtrim(erase(string(X.Substrate),"'"));
    target = double(T.(ecGEM_bioRxns{s}));
    valid = isfinite(target) & target > 0;
    T = T(valid,:);
    target = target(valid);
    modelFields = fieldnames(M);
    gauksFields = modelFields(endsWith(modelFields,'_GAUKS'));
    % Add 13C exchange substrates to the common substrate list.
    if s <= 2
        c13ExchangeRxns = string(C13.RxnName(c13ConstraintRows));
        testG = M.(gauksFields{1});
        modelRxns = string(testG.rxns);
        c13ExchangeRxns = c13ExchangeRxns( ...
            ismember(c13ExchangeRxns+"_basal",modelRxns) & ...
            ismember(c13ExchangeRxns+"_extended",modelRxns));
    else
        c13ExchangeRxns = strings(0,1);
    end
    allCarbon = unique([T.Substrate;X.Substrate;c13ExchangeRxns],'stable');
    allCarbon = allCarbon(~ismissing(allCarbon) & strlength(allCarbon)>0);
    muData = double(X.(ecGEM_bioRxns{s}));
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
        mB.c(:) = 0;
        mG.c(:) = 0;
        mB.c(ibB) = 1;
        mG.c(ibG) = 1;
        idsB = findRxnIDs(mB,cellstr(allCarbon));
        idsG = findRxnIDs(mG,cellstr(allCarbon));
        params = extractSluiceParams(mG,cellstr(allCarbon));
        pool = find(strcmp(mG.mets,'prot_pool'));
        extAll = find(endsWith(string(mG.rxns),"_extended"));
        statusByRow = repmat("NoReference",height(T),1);
        %% Maximum growth, GAUKS activation and parameter audit
        for i = 1:height(T)
            ex = T.Substrate(i);
            ox = string(T.OxAvail(i));
            k = find(allCarbon==ex,1);
            if isempty(k)
                continue
            end
            jB = idsB(k);
            jG = idsG(k);
            iBasal = findRxnIDs(mG,char(ex+"_basal"));
            iExtended = findRxnIDs(mG,char(ex+"_extended"));
            if jB == 0 || jG == 0 || iBasal == 0 || iExtended == 0
                continue
            end
            umin = params.umin(k);
            xi = params.xi(k);
            b = mB;
            g = mG;
            b.lb(idsB(idsB>0)) = 0;
            b.ub(idsB(idsB>0)) = max(0,b.ub(idsB(idsB>0)));
            b.lb(jB) = -uptakeCap;
            b.ub(jB) = 0;
            g.lb(idsG(idsG>0)) = 0;
            g.ub(idsG(idsG>0)) = max(0,g.ub(idsG(idsG>0)));
            g.lb(jG) = -uptakeCap;
            g.ub(jG) = 0;
            if ox == "anaerobic"
                ioB = findRxnIDs(b,oxygenRxns{s});
                ioG = findRxnIDs(g,oxygenRxns{s});
                b = anaerobicModel(b, orgNames{s}, ioB);
                g = anaerobicModel(g, orgNames{s}, ioG);
            end
            sb = solveLP(b,'max');
            sg = solveLP(g,'max');
            muB = NaN;
            muG = NaN;
            qB = NaN;
            qG = NaN;
            if sb.stat == 1
                muB = sb.x(ibB);
                qB = max(0,-sb.x(jB));
            end
            if sg.stat == 1
                muG = sg.x(ibG);
                qG = max(0,-sg.x(jG));
            end
            searchTriggered = isfinite(muB) && muB > target(i)+growthTol;
            % Umin = 1000 means GAUKS was not triggered for this substrate.
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
            elseif ~searchTriggered
                status = "Inactive_NoOverprediction";
            else
                status = "Active";
            end
            active = status == "Active";
            statusByRow(i) = status;
            growthChange = muG-muB;
            growthChangePct = NaN;
            if isfinite(muB) && abs(muB)>growthTol
                growthChangePct = 100*growthChange/muB;
            end
            omegaB = abs(muB-target(i))/target(i);
            omegaG = abs(muG-target(i))/target(i);
            R(end+1,:) = {speciesNames(s),string(baseName),string(gauksName),...
                ex,ox,target(i),muB,muG,growthChange,growthChangePct,...
                qB,qG,qG-qB,omegaB,omegaG,omegaG-omegaB,status};
            P(end+1,:) = {speciesNames(s),string(baseName),string(gauksName),...
                ex,ox,umin,xi,target(i),muB,muG,qB,qG,...
                searchTriggered,status,active};
            %% A: GAUKS threshold scan
            %% A: GAUKS threshold scan
            if active
                for r = ratios
                    u = r*umin;
                    if u > uptakeCap || u <= 0
                        continue
                    end
                    ba = b;
                    ga = g;
                    ba.lb(jB) = -u;
                    ba.ub(jB) = 0;
                    ga.lb(jG) = -u;
                    ga.ub(jG) = 0;
                    saB = solveLP(ba,'max');
                    saG = solveLP(ga,'max');
                    yB = NaN;
                    yG = NaN;
                    basalFlux = NaN;
                    extendedFlux = NaN;
                    proteinExtra = NaN;
                    balanceError = NaN;
                    if saB.stat == 1
                        yB = saB.x(ibB);
                    end
                    if saG.stat == 1
                        yG = saG.x(ibG);
                        basalFlux = max(0,-saG.x(iBasal));
                        extendedFlux = max(0,-saG.x(iExtended));
                        proteinExtra = xi*extendedFlux;
                        balanceError = -(saG.x(iBasal)+saG.x(iExtended))-u;
                    end
                    A(end+1,:) = {speciesNames(s),string(baseName),ex,ox,...
                        umin,xi,r,u,saB.stat==1,saG.stat==1,...
                        yB,yG,yG-yB,basalFlux,extendedFlux,...
                        proteinExtra,balanceError,status};
                end
            end
            %% C2: fixed-growth uptake FVA
            if isfinite(muB) && isfinite(muG) && min(muB,muG)>growthTol
                for f = growthFractions
                    fixedGrowth = f*min(muB,muG);
                    tb = b;
                    tg = g;
                    tb.lb(ibB) = fixedGrowth;
                    tb.ub(ibB) = fixedGrowth;
                    tg.lb(ibG) = fixedGrowth;
                    tg.ub(ibG) = fixedGrowth;
                    tb.c(:) = 0;
                    tb.c(jB) = 1;
                    slB = solveLP(tb,'max');
                    tb.c(jB) = -1;
                    shB = solveLP(tb,'max');
                    tg.c(:) = 0;
                    tg.c(jG) = 1;
                    slG = solveLP(tg,'max');
                    tg.c(jG) = -1;
                    shG = solveLP(tg,'max');
                    uv = NaN(2,2);
                    ok = [false false];
                    ok(1) = slB.stat==1 && shB.stat==1;
                    ok(2) = slG.stat==1 && shG.stat==1;
                    if ok(1)
                        uv(1,:) = [-slB.x(jB),-shB.x(jB)];
                    end
                    if ok(2)
                        uv(2,:) = [-slG.x(jG),-shG.x(jG)];
                    end
                    C(end+1,:) = {speciesNames(s),string(baseName),"C2",NaN,...
                        ex,ox,xi,umin,fixedGrowth,ok(1),ok(2),...
                        uv(1,1),uv(1,2),uv(2,1),uv(2,2),...
                        uv(1,2)-uv(1,1),uv(2,2)-uv(2,1),...
                        diff(uv(2,:))-diff(uv(1,:)),...
                        uv(1,2)>=999.999,uv(2,2)>=999.999,status};
                end
            end
            %% Z: zero-protein-penalty control
            if ~isempty(pool)
                g0 = g;
                g0.S(pool,extAll) = 0;
                for r = [0.5 1 1.5 2]
                    u = r*umin;
                    if u > uptakeCap
                        continue
                    end
                    bb0 = b;
                    gg0 = g0;
                    bb0.lb(jB) = -u;
                    bb0.ub(jB) = -u;
                    gg0.lb(jG) = -u;
                    gg0.ub(jG) = -u;
                    szB = solveLP(bb0,'max');
                    szG = solveLP(gg0,'max');
                    mzB = NaN;
                    mzG = NaN;
                    if szB.stat == 1, mzB = szB.x(ibB); end
                    if szG.stat == 1, mzG = szG.x(ibG); end
                    Z(end+1,:) = {speciesNames(s),string(baseName),ex,ox,...
                        r,mzB,mzG,mzG-mzB,status};
                end
            end
        end
        %% B, C1 and D: Bayesian experimental conditions
        for kx = 1:height(X)
            ox = string(X.OxAvail(kx));
            if ox == "limited"
                continue
            end
            ex = X.Substrate(kx);
            ix = find(allCarbon==ex,1);
            if isempty(ix)
                continue
            end
            jB = idsB(ix);
            jG = idsG(ix);
            iExtended = findRxnIDs(mG,char(ex+"_extended"));
            if jB==0 || jG==0 || iExtended==0
                continue
            end
            tMatch = find(T.Substrate==ex & string(T.OxAvail)==ox,1);
            if isempty(tMatch)
                tMatch = find(T.Substrate==ex,1);
            end
            if isempty(tMatch)
                statusX = "NoReference";
            else
                statusX = statusByRow(tMatch);
            end
            umin = params.umin(ix);
            xi = params.xi(ix);
            u = abs(double(X.Uptake(kx)));
            muExp = muData(kx);
            if ~isfinite(u) || ~isfinite(muExp) || muExp<=0
                continue
            end
            ru = NaN;
            if isfinite(umin) && umin>paramTol
                ru = u/umin;
            end
            grp = "B0";
            if statusX=="Active" && isfinite(ru)
                if ru<=1+paramTol
                    grp = "B1";
                else
                    grp = "B2";
                end
            end
            bX = mB;
            gX = mG;
            bX.lb(idsB(idsB>0)) = 0;
            bX.ub(idsB(idsB>0)) = max(0,bX.ub(idsB(idsB>0)));
            gX.lb(idsG(idsG>0)) = 0;
            gX.ub(idsG(idsG>0)) = max(0,gX.ub(idsG(idsG>0)));
            bX.lb(jB) = -uptakeCap;
            bX.ub(jB) = 0;
            gX.lb(jG) = -uptakeCap;
            gX.ub(jG) = 0;
            if ox=="anaerobic"
                ioB = findRxnIDs(bX,oxygenRxns{s});
                ioG = findRxnIDs(gX,oxygenRxns{s});
                if ioB>0, bX.lb(ioB)=0; end
                if ioG>0, gX.lb(ioG)=0; end
            end
            %% B: fixed experimental uptake, predict growth
            bb = bX;
            gg = gX;
            bb.lb(jB) = -u;
            bb.ub(jB) = -u;
            gg.lb(jG) = -u;
            gg.ub(jG) = -u;
            sb = solveLP(bb,'max');
            sg = solveLP(gg,'max');
            yB = NaN;
            yG = NaN;
            ext = NaN;
            extra = NaN;
            if sb.stat==1
                yB = sb.x(ibB);
            end
            if sg.stat==1
                yG = sg.x(ibG);
                ext = max(0,-sg.x(iExtended));
                extra = xi*ext;
            end
            wB = abs(yB-muExp)/muExp;
            wG = abs(yG-muExp)/muExp;
            B(end+1,:) = {speciesNames(s),string(baseName),kx,ex,ox,...
                umin,xi,ru,grp,u,muExp,sb.stat==1,sg.stat==1,...
                yB,yG,wB,wG,wG-wB,ext,extra,statusX};
            %% C1: fixed experimental growth, optimize uptake range
            tb = bX;
            tg = gX;
            tb.lb(ibB) = muExp;
            tb.ub(ibB) = muExp;
            tg.lb(ibG) = muExp;
            tg.ub(ibG) = muExp;
            tb.c(:) = 0;
            tb.c(jB) = 1;
            slB = solveLP(tb,'max');
            tb.c(jB) = -1;
            shB = solveLP(tb,'max');
            tg.c(:) = 0;
            tg.c(jG) = 1;
            slG = solveLP(tg,'max');
            tg.c(jG) = -1;
            shG = solveLP(tg,'max');
            uv = NaN(2,2);
            ok = [false false];
            ok(1) = slB.stat==1 && shB.stat==1;
            ok(2) = slG.stat==1 && shG.stat==1;
            if ok(1)
                uv(1,:) = [-slB.x(jB),-shB.x(jB)];
            end
            if ok(2)
                uv(2,:) = [-slG.x(jG),-shG.x(jG)];
            end
            C(end+1,:) = {speciesNames(s),string(baseName),"C1",kx,...
                ex,ox,xi,umin,muExp,ok(1),ok(2),...
                uv(1,1),uv(1,2),uv(2,1),uv(2,2),...
                uv(1,2)-uv(1,1),uv(2,2)-uv(2,1),...
                diff(uv(2,:))-diff(uv(1,:)),...
                uv(1,2)>=999.999,uv(2,2)>=999.999,statusX};
            %% D: fixed experimental uptake and growth
            sdB = bb;
            sdG = gg;
            sdB.lb(ibB) = muExp;
            sdB.ub(ibB) = muExp;
            sdG.lb(ibG) = muExp;
            sdG.ub(ibG) = muExp;
            sdBsol = solveLP(sdB,'max');
            sdGsol = solveLP(sdG,'max');
            extD = NaN;
            extraD = NaN;
            if sdGsol.stat==1
                extD = max(0,-sdGsol.x(iExtended));
                extraD = xi*extD;
            end
            D(end+1,:) = {speciesNames(s),string(baseName),kx,ex,ox,...
                u,muExp,ru,sdBsol.stat==1,sdGsol.stat==1,...
                extD,extraD,statusX};
        end
        %% E: 13C internal flux RMSE
        if s <= 2
            c13Status = "NoReference";
            idx13 = find(T.Substrate==cSources{s} & string(T.OxAvail)~="anaerobic" & string(T.OxAvail)~="limited",1);
            if isempty(idx13)
                idx13 = find(T.Substrate==cSources{s},1);
            end
            if ~isempty(idx13)
                c13Status = statusByRow(idx13);
            end
            for cc = 1:numel(c13CondCols)
                condName = string(c13ConditionNames{cc});
                condCol = c13CondCols(cc);
                bC = mB;
                gC = mG;
                bC.c(:) = 0;
                gC.c(:) = 0;
                bC.lb(idsB(idsB>0)) = 0;
                bC.ub(idsB(idsB>0)) = max(0,bC.ub(idsB(idsB>0)));
                gC.lb(idsG(idsG>0)) = 0;
                gC.ub(idsG(idsG>0)) = max(0,gC.ub(idsG(idsG>0)));
                % Apply measured biomass and reaction constraints
                for ci = c13ConstraintRows(:)'
                    rx = strtrim(string(C13.RxnName(ci)));
                    val = double(C13{ci,condCol});
                    if ~isfinite(val)
                        continue
                    end
                    if rx == string(ecGEM_bioRxns{s})
                        idB = findRxnIDs(bC,ecGEM_bioRxns{s});
                        idG = findRxnIDs(gC,ecGEM_bioRxns{s});
                        bC.lb(idB) = val;
                        bC.ub(idB) = val;
                        gC.lb(idG) = val;
                        gC.ub(idG) = val;
                    else
                        idB = findRxnIDs(bC,char(rx));
                        idG = findRxnIDs(gC,char(rx));
                        if idB>0
                            bC.lb(idB) = val;
                            bC.ub(idB) = val;
                        end
                        if idG>0
                            gC.lb(idG) = val;
                            gC.ub(idG) = val;
                        end
                    end
                end
                % Select representative flux solutions using the protein-pool objective
                pB = findRxnIDs(bC,'prot_pool_exchange');
                pG = findRxnIDs(gC,'prot_pool_exchange');
                bC.c(:) = 0;
                gC.c(:) = 0;
                bC.c(pB) = 1;
                gC.c(pG) = 1;
                sbC = solveLP(bC);
                sgC = solveLP(gC);
                feasibleB = sbC.stat==1 && ~isempty(sbC.x);
                feasibleG = sgC.stat==1 && ~isempty(sgC.x);
                expFlux = [];
                simFluxB = [];
                simFluxG = [];
                carbonNums = [];
                % Compare measured and simulated internal fluxes
                if feasibleB && feasibleG
                    for fr = c13FluxRows(:)'
                        rxnRaw = strtrim(string(C13.RxnName(fr)));
                        dirRaw = strtrim(string(C13.Direction(fr)));
                        vExp = double(C13{fr,condCol});
                        if ~isfinite(vExp)
                            continue
                        end
                        if startsWith(rxnRaw,"{")
                            rxnList = strsplit(char(erase(erase(rxnRaw,"{"),"}")),';');
                        else
                            rxnList = {char(rxnRaw)};
                        end
                        rxnList = strtrim(rxnList);
                        dirStr = char(erase(erase(dirRaw,"{"),"}"));
                        dirVals = str2double(strsplit(dirStr,';'));
                        if isempty(dirVals) || any(isnan(dirVals))
                            dirVals = ones(1,numel(rxnList));
                        end
                        if numel(dirVals)==1 && numel(rxnList)>1
                            dirVals = repmat(dirVals,1,numel(rxnList));
                        end
                        dirVals = dirVals(1:numel(rxnList));
                        dirVals = dirVals(:);
                        idsFB = findRxnIDs(bC,rxnList);
                        idsFG = findRxnIDs(gC,rxnList);
                        idsFB = idsFB(:);
                        idsFG = idsFG(:);
                        if numel(idsFB)~=numel(rxnList) || numel(idsFG)~=numel(rxnList) || any(idsFB==0) || any(idsFG==0)
                            continue
                        end
                        expFlux(end+1,1) = vExp;
                        simFluxB(end+1,1) = sum(dirVals.*sbC.x(idsFB));
                        simFluxG(end+1,1) = sum(dirVals.*sgC.x(idsFG));
                        carbonNums(end+1,1) = str2double(string(C13.Carbon(fr)));
                    end
                end
                rmseFluxB = NaN;
                rmseFluxG = NaN;
                rmseCarbonB = NaN;
                rmseCarbonG = NaN;
                if ~isempty(expFlux)
                    rmseFluxB = sqrt(mean((expFlux-simFluxB).^2));
                    rmseFluxG = sqrt(mean((expFlux-simFluxG).^2));
                    validCarbon = isfinite(carbonNums);
                    if any(validCarbon)
                        rmseCarbonB = sqrt(mean(((expFlux(validCarbon)-simFluxB(validCarbon)).*carbonNums(validCarbon)).^2));
                        rmseCarbonG = sqrt(mean(((expFlux(validCarbon)-simFluxG(validCarbon)).*carbonNums(validCarbon)).^2));
                    end
                end
                E(end+1,:) = {speciesNames(s),string(baseName),string(gauksName),...
                    condName,cSources{s},rmseFluxB,rmseFluxG,rmseFluxG-rmseFluxB,...
                    rmseCarbonB,rmseCarbonG,rmseCarbonG-rmseCarbonB,...
                    numel(expFlux),sum(isfinite(carbonNums)),...
                    feasibleB,feasibleG,c13Status};
            end
        end
    end
end
%% Save all results into one Excel workbook
outFile = fullfile(findECOMAProot, 'projects', ...
    'GAUKS_Combined_Results.xlsx');
if exist(outFile, 'file')
    delete(outFile);
end
% Convert all results into tables
resultTables = {
    cell2table(R, 'VariableNames', { ...
        'Species','BaseModel','GAUKSModel','Substrate','OxAvail', ...
        'GrowthReference','GrowthBefore','GrowthAfter','GrowthChange', ...
        'GrowthChangePercent','UptakeBefore','UptakeAfter','UptakeChange', ...
        'OmegaBefore','OmegaAfter','DeltaOmega','GAUKSStatus'})
    cell2table(P, 'VariableNames', { ...
        'Species','BaseModel','GAUKSModel','Substrate','OxAvail', ...
        'Umin','Xi','MuTarget','MuBefore','MuAfter','UptakeBefore', ...
        'UptakeAfter','SearchTriggered','GAUKSStatus','GAUKSActive'})
    cell2table(A, 'VariableNames', { ...
        'Species','Model','Substrate','OxAvail','Umin','Xi','RU', ...
        'UptakeCap','FeasibleBefore','FeasibleAfter','MuBefore', ...
        'MuAfter','DeltaMu','Basal','Extended','ProteinExtra', ...
        'BalanceError','GAUKSStatus'})
    cell2table(B, 'VariableNames', { ...
        'Species','Model','RecordID','Substrate','OxAvail','Umin','Xi','RU', ...
        'Group','UExp','MuExp','FeasibleBefore','FeasibleAfter','MuBefore', ...
        'MuAfter','OmegaBefore','OmegaAfter','DeltaOmega','Extended', ...
        'ProteinExtra','GAUKSStatus'})
    cell2table(C, 'VariableNames', { ...
        'Species','Model','Test','RecordID','Substrate','OxAvail','Xi','Umin', ...
        'GrowthFixed','FeasibleBefore','FeasibleAfter','UminBefore', ...
        'UmaxBefore','UminAfter','UmaxAfter','DeltaUmin','DeltaUmax', ...
        'DeltaWidth','BeforeHitCap','AfterHitCap','GAUKSStatus'})
    cell2table(D, 'VariableNames', { ...
        'Species','Model','RecordID','Substrate','OxAvail','UExp','MuExp', ...
        'RU','FeasibleBefore','FeasibleAfter','Extended','ProteinExtra', ...
        'GAUKSStatus'})
    cell2table(E, 'VariableNames', { ...
        'Species','BaseModel','GAUKSModel','Condition','CarbonSource', ...
        'RMSE_Flux_Before','RMSE_Flux_After','DeltaRMSE_Flux', ...
        'RMSE_Carbon_Before','RMSE_Carbon_After','DeltaRMSE_Carbon', ...
        'NMatchedFlux','NMatchedCarbon','FeasibleBefore','FeasibleAfter', ...
        'GAUKSStatus'})
    cell2table(Z, 'VariableNames', { ...
        'Species','Model','Substrate','OxAvail','RU','MuBefore', ...
        'MuXiZero','DeltaMu','GAUKSStatus'})
};
sheetNames = { ...
    'MaxGrowth', ...
    'ParamAudit', ...
    'ThresholdScan', ...
    'ExperimentalGrowth', ...
    'UptakeFVA', ...
    'PhysFeasibility', ...
    'Internal13CRMSE', ...
    'ZeroPenalty'};
% Write all result tables into the same workbook
for k = 1:numel(resultTables)
    writetable(resultTables{k}, outFile, ...
        'Sheet', sheetNames{k});
end
fprintf('All GAUKS analyses saved to:\n%s\n', outFile);



