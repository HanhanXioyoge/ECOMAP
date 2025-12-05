%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,P,C,R,D,L] = sumBioMass(model)
% Calculates breakdown of biomass for the yeast model:
% X -> Biomass fraction without lipids [g/gDW]
% P -> Protein fraction [g/gDW]
% C -> Carbohydrate fraction [g/gDW]
% R -> RNA fraction [g/gDW]
% D -> DNA fraction [g/gDW]
% L -> Lipid fraction [g/gDW]
% Function adapted from SLIMEr: https://github.com/SysBioChalmers/SLIMEr
%
% Benjamin Sanchez. Last update: 2018-10-23
% Ivan Domenzain.   Last update: 2019-07-13
% 25.05.2025 Qilin Chen, Adjusted to fit ECOMAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,P,C,R,D,L] = sumBioMass(model, parameters)

if nargin < 5 || isempty(parameters)
    parameters = ParameterManager.getParams();
    if isempty(parameters), error('ParameterManager is not set.'); end
end

comps     = parameters.PRESTO.comps;
bioRxn    = parameters.bioRxn;

%Get main fractions:
[P,X] = getFraction(model,comps,'P',0);
[C,X] = getFraction(model,comps,'C',X);
[R,X] = getFraction(model,comps,'R',X);
[D,X] = getFraction(model,comps,'D',X);
[L,X] = getFraction(model,comps,'L',X);

%Add up any remaining components:
bioPos = strcmp(model.rxns,bioRxn);
for i = 1:length(model.mets)
    pos = strcmp(comps(:,1),model.mets{i});
    if sum(pos) == 1
        abundance = -model.S(i,bioPos)*comps{pos,2}/1000;
        X         = X + abundance;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,X] = getFraction(model,comps,compType,X)

%Define pseudoreaction name:
rxnName = [compType ' pseudoreaction'];
rxnName = strrep(rxnName,'P','protein');
rxnName = strrep(rxnName,'C','carbohydrate');
rxnName = strrep(rxnName,'R','RNA');
rxnName = strrep(rxnName,'D','DNA');
rxnName = strrep(rxnName,'L','lipid backbone');

%Add up fraction:
fractionPos = strcmp(model.rxnNames,rxnName);
if contains(rxnName,'lipid')
    subs = model.S(:,fractionPos) < 0;        %substrates in pseudo-rxn
    F    = -sum(model.S(subs,fractionPos));   %g/gDW
else
    comps = comps(strcmp(comps(:,3),compType),:);
    F = 0;
    %Add up all components:
    for i = 1:length(model.mets)
        pos = strcmp(comps(:,1),strtok(model.mets{i},'['));
        if sum(pos) == 1
            abundance = -model.S(i,fractionPos)*(comps{pos,2}-18)/1000;
            F         = F + abundance;
        end
    end
end
X = X + F;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%