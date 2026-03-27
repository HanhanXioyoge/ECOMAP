function model = anaerobicModel(model, org_name, EX_o2)

if nargin < 3
    EX_o2 = ''; 
end

if strcmp(org_name, 'Saccharomyces cerevisiae')
    GAM   = 58.1988;
    NGAM  = 0;
    model = changeGAM(model,GAM,NGAM);
    mets = {'s_3714','s_1198','s_1203','s_1207','s_1212','s_0529'};
    [~,met_index] = ismember(mets,model.mets);
    model.S(met_index,strcmp(model.rxns,'r_4598')) = 0;
    
    model.lb(strcmp(model.rxns,'r_1992')) = 0;        %O2
    model.lb(strcmp(model.rxns,'r_1757')) = -1000;    %ergosterol
    model.lb(strcmp(model.rxns,'r_1915')) = -1000;    %lanosterol
    model.lb(strcmp(model.rxns,'r_1994')) = -1000;    %palmitoleate
    model.lb(strcmp(model.rxns,'r_2106')) = -1000;    %zymosterol
    model.lb(strcmp(model.rxns,'r_2134')) = -1000;    %14-demethyllanosterol
    model.lb(strcmp(model.rxns,'r_2137')) = -1000;    %ergosta-5,7,22,24(28)-tetraen-3beta-ol
    model.lb(strcmp(model.rxns,'r_2189')) = -1000;    %oleate
    model = changeRxnBounds(model,'r_2090',-1000,'l');% uracil uptake
    %4rd change: Blocked pathways for proper glycerol production
    %Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
    model.lb(startsWith(model.rxns,'r_0713')) = 0; %Mithocondria
    model.lb(strcmp(model.rxns,'r_0714')) = 0; %Cytoplasm
    model.lb(strcmp(model.rxns,'r_0713_rvs')) = 0; %Mithocondria
    model.lb(strcmp(model.rxns,'r_0714')) = 0; %Cytoplasm
    model.lb(strcmp(model.rxns,'r_0714_rvs')) = 0; %Cytoplasm
    %Block glycerol dehydroginase (only acts in microaerobic conditions)
    model.ub(strcmp(model.rxns,'r_0487')) = 0;
    model.ub(strcmp(model.rxns,'r_0487_rvs')) = 0;
    %Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
    model.ub(strcmp(model.rxns,'r_0472')) = 0;
    model.ub(strcmp(model.rxns,'r_0472_fwd')) = 0;
    
    %4th change: Blocked pathways for proper glycerol production
    %Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
    idx = regexp(cellstr(model.rxns),'r_0713[\w*]rvs');
    idx = find(cellfun(@isempty,idx)==0);
    model = changeRxnBounds(model,model.rxns(idx),0,'b'); 
    model.lb(strcmp(model.rxns,'r_0713')) = 0; %Mithocondria % in case this one does not have any grRule
    
    idx = regexp(cellstr(model.rxns),'r_0714[\w*]rvs');
    idx = find(cellfun(@isempty,idx)==0);
    model = changeRxnBounds(model,model.rxns(idx),0,'b'); 
    model.lb(strcmp(model.rxns,'r_0714')) = 0; %Cytoplasm
    
    %Block glycerol dehydroginase (only acts in microaerobic conditions)
    % model.ub(strcmp(model.rxns,'r_0487')) = 0; is blocked in blockrxns
    
    %Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
    idx = find(startsWith(model.rxns,'r_0472_'));
    model = changeRxnBounds(model,model.rxns(idx),0,'b');
    model.ub(strcmp(model.rxns,'r_0472')) = 0;
elseif strcmp(org_name, 'Escherichia coli')
    model.lb(strcmp(model.rxns, 'EX_o2_e')) = 0;
elseif strcmp(org_name, 'Corynebacterium glutamicum')
    model.lb(strcmp(model.rxns, 'EX_o2_e')) = 0;
else
    model.lb(strcmp(model.rxns, EX_o2)) = 0;
end
end

function model = changeGAM(model,GAM,NGAM)

bioPos = strcmp(model.rxnNames,'biomass pseudoreaction');
for i = 1:length(model.mets)
    S_ix  = model.S(i,bioPos);
    isGAM = sum(strcmp({'ATP [cytoplasm]','ADP [cytoplasm]','H2O [cytoplasm]', ...
        'H+ [cytoplasm]','phosphate [cytoplasm]'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAM
        model.S(i,bioPos) = sign(S_ix)*GAM;
    end
end

if nargin >1
    pos = strcmp(model.rxnNames,'non-growth associated maintenance reaction');%NGAM
    %model = setParam(model,'eq',model.rxns(pos),NGAM);% set both lb and ub to be mu
    model.lb(pos) = NGAM;
    model.ub(pos) = NGAM;
    
end

end