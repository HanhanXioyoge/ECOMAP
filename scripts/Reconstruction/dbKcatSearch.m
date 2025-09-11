function kcatList = dbKcatSearch(model, parameters, forceWClvl)
% dbKcatSearch
% Retrieve per-reaction kcat values for an enzyme-constrained model by
% fuzzy-matching EC numbers and substrates against the BRENDA database,
% with fallbacks across substrate specificity, organism proximity, and
% EC-wildcard levels.
%
% SYNTAX
%   kcatList = dbKcatSearch(model, parameters, forceWClvl)
%
% INPUTS
%   model      (struct)  A COBRA/GECKO ecModel with:
%                        - model.enzymeConstraints.rxns      : ec reaction IDs
%                        - model.enzymeConstraints.eccodes   : EC number(s) per ec rxn
%                        - model.enzymeConstraints.ecModeltype (optional): 'basic' or others
%                        Also uses model.S, model.rxns, model.metNames.
%
%   parameters (struct)  Optional. If omitted/empty, uses
%                        ParameterManager.getParams(). Must contain:
%                        - org_name : scientific name used for organism
%                                     matching (e.g., 'Escherichia coli').
%
%   forceWClvl (scalar)  Optional. Default = 0. If >0, progressively relaxes
%                        EC numbers by replacing trailing levels with '-' before
%                        querying (coarser wildcard match). Use sparingly.
%
% OUTPUT
%   kcatList     (struct)  Selected kcat values and provenance:
%       .source       : 'brenda'
%       .rxns         : ec reactions considered (subset of enzyme-constrained rxns)
%       .substrates   : cell array; for each rxn, the substrate names parsed from S
%       .kcats        : numeric vector of chosen kcat [s^-1] (0 if not found)
%       .eccodes      : EC numbers used in the search (after any wildcarding)
%       .wildcardLvl  : per-rxn integer indicating how far the EC had to be relaxed
%       .origin       : per-rxn integer provenance code:
%                       1 = org_s    (same organism, substrate-specific)
%                       2 = rest_s   (other organism, substrate-specific)
%                       3 = org_ns   (same organism, non-substrate-specific)
%                       4 = rest_ns  (other organism, non-substrate-specific)
%                       5 = org_sa   (same organism, from specific activity)
%                       6 = rest_sa  (other organism, from specific activity)
%
% METHOD (high-level)
%   1) Determine the enzyme-constrained reaction set and extract, for each,
%      the substrates and stoichiometric coefficients from model.S.
%   2) Load BRENDA tables of kcat and specific activity (loadBRENDA) and a KEGG-based
%      taxonomy structure (KEGG_struct, find_inKEGG) to estimate phylogenetic proximity.
%   3) Optionally relax EC specificity by adding wildcards (forceWClvl).
%   4) For each reaction, call iterativeMatch to:
%        - prioritize substrate-specific kcats from the same organism,
%        - fall back to close taxa, then to non-substrate-specific values,
%        - finally allow conversion from specific activity (SA) when needed,
%        - record the wildcard level and origin category chosen.
%      An EC→row index is precomputed to accelerate repeated lookups.
%   5) Aggregate per-rxn results into kcatList.
%
% ASSUMPTIONS & NOTES
%   - Units: kcat in s^-1. SA values (e.g., U/mg) are converted internally when used.
%   - If model.enzymeConstraints.ecModeltype == 'basic', the function trims any
%     leading 'ecs_'-style prefix from ec reaction IDs before mapping back to model.rxns.
%   - Phylogenetic distances are approximated from KEGG taxonomy; matches from the
%     target organism are preferred over related organisms.
%   - External utilities required: loadBRENDA, KEGG_struct, find_inKEGG,
%     iterativeMatch, and progressbar.
%
% ERRORS
%   - Missing ParameterManager or parameters.org_name.
%   - Missing model.enzymeConstraints.eccodes.
%
% EXAMPLE
%   params = ParameterManager.getParams();
%   params.org_name = 'Escherichia coli';
%   kcatList = dbKcatSearch(ecModel, params, 0);
%
% SEE ALSO
%   loadBRENDA, iterativeMatch, KEGG_struct, find_inKEGG

    if nargin < 3 || isempty(forceWClvl)
        forceWClvl = 0;
    end
    
    if nargin < 2 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.');
        end
    end

    if ~isfield(model,'enzymeConstraints') || ~isfield(model.enzymeConstraints,'eccodes')
        error('No EC codes in model.enzymeConstraints.eccodes.');
    end

    if ~isfield(parameters,'org_name') || isempty(parameters.org_name)
        error('parameters.org_name is required (e.g., ''Escherichia coli'').');
    end


    if ~isfield(model.enzymeConstraints, 'ecModeltype') || ~strcmp(model.enzymeConstraints.ecModeltype, 'basic')
        rxnNames = model.enzymeConstraints.rxns;
    else
        rxnNames = extractAfter(model.enzymeConstraints.rxns,4);
    end

    ecRxns       = true(numel(model.enzymeConstraints.rxns),1);
    ecRxns       = find(ecRxns);
    org_name     = parameters.org_name;
    eccodes      = model.enzymeConstraints.eccodes(ecRxns);
    substrates   = cell(numel(ecRxns),1);
    substrCoeffs = cell(numel(ecRxns),1);

    [~,originalRxns] = ismember(rxnNames(ecRxns),model.rxns);
    for i = 1:length(ecRxns)
        sel = find(model.S(:,originalRxns(i)) < 0);
        substrates{i}  = model.metNames(sel); 
        substrCoeffs{i} = -model.S(sel,originalRxns(i));
    end

    [KCATcell, SAcell] = loadBRENDA();

    phylDistStruct =  KEGG_struct();

    org_index      = find_inKEGG(org_name,phylDistStruct.names);

    phylDistStruct.genus = lower(regexprep(phylDistStruct.names,'\s.*',''));
    %create a map for the genuses
    phylDistStruct.uniqueGenusList = unique(phylDistStruct.genus);
    phylDistStruct.genusHashMap = containers.Map(phylDistStruct.uniqueGenusList,1:length(phylDistStruct.uniqueGenusList));
    phylDistStruct.uniqueGenusIndices = cell(length(phylDistStruct.uniqueGenusList),1);

    for i = 1:length(phylDistStruct.genus)
        matchInd = cell2mat(values(phylDistStruct.genusHashMap,phylDistStruct.genus(i)));
        phylDistStruct.uniqueGenusIndices{matchInd} = [phylDistStruct.uniqueGenusIndices{matchInd};i];
    end

    %Allocate output
    kcats = zeros(length(eccodes),1);
    mM = length(eccodes);

    %Create empty kcatInfo
    %Legacy, no longer given as output, rather used to construct
    %kcatList.wildcardLvl and kcatList.origin.
    kcatInfo.info.org_s   = zeros(mM,1);
    kcatInfo.info.rest_s  = zeros(mM,1);
    kcatInfo.info.org_ns  = zeros(mM,1);
    kcatInfo.info.rest_ns = zeros(mM,1);
    kcatInfo.info.org_sa  = zeros(mM,1);
    kcatInfo.info.rest_sa = zeros(mM,1);
    kcatInfo.info.wcLevel = NaN(mM,1);
    kcatInfo.stats.queries  = 0;
    kcatInfo.stats.org_s    = 0;
    kcatInfo.stats.rest_s   = 0;
    kcatInfo.stats.org_ns   = 0;
    kcatInfo.stats.rest_ns  = 0;
    kcatInfo.stats.org_sa   = 0;
    kcatInfo.stats.rest_sa  = 0;
    kcatInfo.stats.wc0      = 0;
    kcatInfo.stats.wc1      = 0;
    kcatInfo.stats.wc2      = 0;
    kcatInfo.stats.wc3      = 0;
    kcatInfo.stats.wc4      = 0;
    kcatInfo.stats.matrix   = zeros(6,5);

    %build an EC index to speed things up a bit - many of the ECs appear
    %many times - unnecessary to compare them all
    %so, here, each EC string appears only once, and you get a vector with
    %indices to the rows in KCATcell
    [ECIndexIds,~,ic] = unique(KCATcell{1});
    EcIndexIndices = cell(length(ECIndexIds),1);
    for i = 1:length(EcIndexIndices)
        EcIndexIndices{i} = find(ic == i).';
    end

    %Apply force wildcard level
    while forceWClvl > 0
        eccodes=regexprep(eccodes,'(.)*(\.\d+)(\.-)*$','$1\.-$3');
        forceWClvl = forceWClvl - 1;
    end
    if forceWClvl == 1
        eccodes = regexprep(eccodes,'.*','-\.-\.-\.-');
    end

    progressbar('Gathering kcat values by fuzzy matching to BRENDA database')

    %Main loop:
    for i = 1:mM
        %Match:
        EC = eccodes{i};
        if ~isempty(EC)
            EC = strsplit(EC,';');
            %Try to match direct reaction:
            if ~isempty(substrates{i})
                [kcats(i), kcatInfo.info,kcatInfo.stats] = iterativeMatch(EC,substrates{i},substrCoeffs{i},i,KCATcell,...
                    kcatInfo.info,kcatInfo.stats,org_name,...
                    phylDistStruct,org_index,SAcell,ECIndexIds,EcIndexIndices);
            end
        end
        progressbar(i/mM)
    end

    kcatList.source      = 'brenda';
    kcatList.rxns        = model.enzymeConstraints.rxns(ecRxns);
    kcatList.substrates  = substrates;
    kcatList.kcats       = kcats;
    kcatList.eccodes     = eccodes;
    kcatList.wildcardLvl = kcatInfo.info.wcLevel;
    kcatList.origin      = NaN(numel(model.enzymeConstraints.rxns(ecRxns)),1);

    % This can be refactored, iterativeMatch and their nested functions can
    % just directly report the origin number.
    origin = [kcatInfo.info.org_s kcatInfo.info.rest_s kcatInfo.info.org_ns kcatInfo.info.rest_ns kcatInfo.info.org_sa kcatInfo.info.rest_sa];
    for i=1:6
        kcatList.origin(find(origin(:,i))) = i;
    end
end

%%
function [KCATcell, SAcell] = loadBRENDA()

basePath       = fullfile(findECOMAProot, 'scripts\database');
KCAT_file      = fullfile(basePath,'max_KCAT.txt');
SA_file        = fullfile(basePath,'max_SA.txt');
MW_file        = fullfile(basePath,'max_MW.txt');

%Extract BRENDA DATA from files information
KCATcell       = openDataFile(KCAT_file,1);
scalingFactor = 1/60;    %[umol/min/mg] -> [mmol/s/g]    Old: 60 [umol/min/mg] -> [mmol/h/g]
SA            = openDataFile(SA_file,scalingFactor);
scalingFactor = 1/1000;  %[g/mol] -> [g/mmol]
MW            = openDataFile(MW_file,scalingFactor);

for i=1:4
    SAcell{i} = [];
end
previousEC = []; EC_indexes = [];

%build an index on MW{1} to speed things up a bit
%first just extract the genus (i.e. the first part of the name)
MWECNum = upper(unique(MW{1}));
MWECNumIndices = cell(length(MWECNum),1);
MWECNumHashMap = containers.Map(MWECNum,1:length(MWECNum));
for i = 1:length(MW{1})
    matchInd = cell2mat(values(MWECNumHashMap, MW{1}(i)));
    MWECNumIndices{matchInd} = [MWECNumIndices{matchInd};i];
end


for i=1:length(SA{1})
    %Gets the indexes of the EC repetitions in the MW cell for every
    %new (different) EC
    if ~strcmpi(SA{1}(i), previousEC)
        key = upper(SA{1}(i));
        if isKey(MWECNumHashMap,key) %annoyingly, this seems to be needed
            matchInd = cell2mat(values(MWECNumHashMap,key));
            EC_indexes = MWECNumIndices{matchInd};
        else
            EC_indexes = [];
        end
    end
    mwEC{1} = MW{3}(EC_indexes); mwEC{2} = MW{4}(EC_indexes);
    % just looks for the first match because just the maximal value for
    % each EC# / Orgaism is reported on the file
    org_index = find(strcmpi(SA{3}(i),mwEC{1}),1);
    if ~isempty(org_index)
        SAcell{1} = [SAcell{1};SA{1}(i)];
        SAcell{2} = [SAcell{2};SA{3}(i)];
        SAcell{3} = [SAcell{3}; SA{4}(i)*mwEC{2}(org_index)]; %[1/hr]
        SAcell{4} = [SAcell{4}; mwEC{2}(org_index)];
    end
    previousEC = SA{1}(i);
end

%remove EC in front of all the EC numbers
if ~isempty(KCATcell{1})
    KCATcell{1} = extractAfter(KCATcell{1},2);
end
if ~isempty(SAcell{1})
    SAcell{1} = extractAfter(SAcell{1},2);
end

function data_cell = openDataFile(fileName,scalingFactor)
fID          = fopen(fileName);
data_cell    = textscan(fID,'%q %q %q %f %q','delimiter','\t');
fclose(fID);
data_cell{4} = data_cell{4}*scalingFactor;
%Split string for each organism in the BRENDA data
%{name, taxonomy, KEGG code}
data_cell{3}  = regexprep(data_cell{3},'\/\/.*','');
end
end

%%
function phylDistStruct =  KEGG_struct()
load(fullfile(findRAVENroot(),'external','kegg','keggPhylDist.mat'));
phylDistStruct.ids   = transpose(phylDistStruct.ids);
phylDistStruct.names = transpose(phylDistStruct.names);
phylDistStruct.names = regexprep(phylDistStruct.names,'\s*\(.*','');
end

%%
function org_index = find_inKEGG(org_name,names)
org_index      = find(strcmpi(org_name,names));
if numel(org_index)>1
    org_index = org_index(1);
elseif isempty(org_index)
    org_name    = regexprep(org_name,'\s.*','');
    org_index   = find(strcmpi(org_name,names));
    if numel(org_index)>1
        org_index = org_index(1);
    elseif isempty(org_index)
        org_index = '*';
    end
end
end

%%
function [kcat,dir,tot] =iterativeMatch(EC,subs,substrCoeff,i,KCATcell,dir,tot,...
    name,phylDist,org_index,SAcell,ECIndexIds,EcIndexIndices)
%Will iteratively try to match the EC number to some registry in BRENDA,
%using each time one additional wildcard.

kcat    = zeros(size(EC));
origin  = zeros(size(EC));
matches = zeros(size(EC));
wc_num  = ones(size(EC)).*1000;
for k = 1:length(EC)
    success  = false;
    while ~success
        %Atempt match:
        [kcat(k),origin(k),matches(k)] = mainMatch(EC{k},subs,substrCoeff,KCATcell,...
            name,phylDist,...
            org_index,SAcell,ECIndexIds,EcIndexIndices);
        %If any match found, ends. If not, introduces one extra wild card and
        %tries again:
        if origin(k) > 0
            success   = true;
            wc_num(k) = sum(EC{k}=='-');
        else
            dot_pos  = [2 strfind(EC{k},'.')];
            wild_num = sum(EC{k}=='-');
            wc_text  = '-.-.-.-';
            EC{k}    = [EC{k}(1:dot_pos(4-wild_num)) wc_text(1:2*wild_num+1)];
        end
    end
end

if sum(origin) > 0
    %For more than one EC: Choose the maximum value among the ones with the
    %less amount of wildcards and the better origin:
    best_pos   = (wc_num == min(wc_num));
    new_origin = origin(best_pos);
    best_pos   = (origin == min(new_origin(new_origin~=0)));
    max_pos    = find(kcat == max(kcat(best_pos)));
    wc_num     = wc_num(max_pos(1));
    origin     = origin(max_pos(1));
    matches    = matches(max_pos(1));
    kcat       = kcat(max_pos(1));

    %Update dir and tot:
    dir.org_s(i)   = matches*(origin == 1);
    dir.rest_s(i)  = matches*(origin == 2);
    dir.org_ns(i)  = matches*(origin == 3);
    dir.org_sa(i)  = matches*(origin == 4);
    dir.rest_ns(i) = matches*(origin == 5);
    dir.rest_sa(i) = matches*(origin == 6);
    dir.wcLevel(i) = wc_num;
    tot.org_s        = tot.org_s   + (origin == 1);
    tot.rest_s       = tot.rest_s  + (origin == 2);
    tot.org_ns       = tot.org_ns  + (origin == 3);
    tot.org_sa       = tot.org_sa  + (origin == 4);
    tot.rest_ns      = tot.rest_ns + (origin == 5);
    tot.rest_sa      = tot.rest_sa + (origin == 6);
    tot.wc0          = tot.wc0     + (wc_num == 0);
    tot.wc1          = tot.wc1     + (wc_num == 1);
    tot.wc2          = tot.wc2     + (wc_num == 2);
    tot.wc3          = tot.wc3     + (wc_num == 3);
    tot.wc4          = tot.wc4     + (wc_num == 4);
    tot.queries      = tot.queries + 1;
    tot.matrix(origin,wc_num+1) = tot.matrix(origin,wc_num+1) + 1;
end

end

%%
function [kcat,origin,matches] = mainMatch(EC,subs,substrCoeff,KCATcell,...
    name,phylDist,org_index,SAcell,ECIndexIds,EcIndexIndices)

%First make the string matching. This takes time, so we only want to do
%this once:
%Relaxes matching if wild cards are present:
wild     = false;
wild_pos = strfind(EC,'-');
if ~isempty(wild_pos)
    EC   = EC(1:wild_pos(1)-1);
    wild = true;
end
stringMatchesEC_cell = extract_string_matches(EC,KCATcell{1},wild,ECIndexIds,EcIndexIndices);

% Matching function prioritizing organism and substrate specificity when
% available.

origin = 0;
%First try to match organism and substrate:
[kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,name,true,false,...
    phylDist,org_index,SAcell,stringMatchesEC_cell,[]);
if matches > 0 && ~wild % If wildcard, ignore substrate match
    origin = 1;
    %If no match, try the closest organism but match the substrate:
else
    [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,'',true,false,...
        phylDist,org_index,SAcell,stringMatchesEC_cell,[]);
    if matches > 0 && ~wild % If wildcard, ignore substrate match
        origin = 2;
        %If no match, try to match organism but with any substrate:
    else
        [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,name,false,false,...
            phylDist,org_index,SAcell,stringMatchesEC_cell,[]);
        if matches > 0
            origin = 3;
            %If no match, try to match organism but for any substrate (SA*MW):
        else
            %create matching index for SA, has not been needed until now
            stringMatchesSA = extract_string_matches(EC,SAcell{1},wild,[],[]);

            [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,name,false,...
                true,phylDist,org_index,...
                SAcell,stringMatchesEC_cell,stringMatchesSA);
            if matches > 0
                origin = 4;
                %If no match, try any organism and any substrate:
            else
                [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,'',false,...
                    false,phylDist,...
                    org_index,SAcell,stringMatchesEC_cell,stringMatchesSA);
                if matches > 0
                    origin = 5;
                    %Again if no match, look for any org and SA*MW
                else
                    [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,'',...
                        false,true,phylDist,...
                        org_index,SAcell,stringMatchesEC_cell,stringMatchesSA);
                    if matches > 0
                        origin = 6;
                    end
                end

            end
        end
    end
end
end

%%
function [kcat,matches] = matchKcat(EC,subs,substrCoeff,KCATcell,organism,...
    substrate,SA,phylDist,...
    org_index,SAcell,KCATcellMatches,SAcellMatches)

%Will go through BRENDA and will record any match. Afterwards, it will
%return the average value and the number of matches attained.
kcat    = [];
matches = 0;

if SA
    %SAcell{1},wild,[],[]
    EC_indexes = extract_indexes(SAcellMatches,[],SAcell{2},subs,substrate,...
        organism,org_index,phylDist);

    kcat       = SAcell{3}(EC_indexes);
    org_cell   = SAcell{2}(EC_indexes);
    MW_BRENDA  = SAcell{4}(EC_indexes);

else
    %KCATcell{1},wild,ECIndexIds,EcIndexIndices
    EC_indexes = extract_indexes(KCATcellMatches,KCATcell{2},KCATcell{3},...
        subs,substrate,organism,org_index,...
        phylDist);
    if substrate
        for j = 1:length(EC_indexes)
            indx = EC_indexes(j);
            for k = 1:length(subs)
                if (isempty(subs{k}))
                    break;
                end
                %l = logical(strcmpi(model.metNames,subs{k}).*(model.S(:,i)~=0)); %I don't understand the .* (model.S(:,i)~=0) part, it shouldn't be needed?/JG;
                if ~isempty(subs{k}) && strcmpi(subs{k},KCATcell{2}(indx))
                    if KCATcell{4}(indx) > 0
                        coeff = min(substrCoeff);
                        kCatTmp = KCATcell{4}(indx);
                        kcat  = [kcat;kCatTmp/coeff];
                    end
                end
            end
        end
    else
        kcat = KCATcell{4}(EC_indexes);
    end
end
%Return maximum value:
if isempty(kcat)
    kcat = 0;
else
    matches        = length(kcat);
    [kcat,MaxIndx] = max(kcat);
end
%Avoid SA*Mw values over the diffusion limit rate  [Bar-Even et al. 2011]
if kcat>(1E7)
    kcat = 1E7;
end
end

%%
%Make the string matches of the ECs. This is heavy, so only do it once!
function EC_indexes = extract_string_matches(EC,EC_cell,wild,ECIndexIds,EcIndexIndices)
EC_indexes = [];
EC_indexesOld = [];
if wild
    if (~isempty(ECIndexIds)) %In some cases the EC_cell is not from KCatCell
        X = find(contains(ECIndexIds, EC));
        for j = 1:length(X)
            EC_indexes = [EC_indexes,EcIndexIndices{X(j)}];
        end
    else %Not optimized
        for j=1:length(EC_cell)
            if strfind(EC_cell{j},EC)==1
                EC_indexes = [EC_indexes,j];
            end
        end
    end
else
    if (~isempty(ECIndexIds)) %In some cases the EC_cell is not from KCatCell
        mtch = find(strcmpi(EC,ECIndexIds));
        if (~isempty(mtch))
            EC_indexes = EcIndexIndices{mtch};
        end
    else %%Not optimized
        if ~isempty(EC_cell)
            EC_indexes = transpose(find(strcmpi(EC,EC_cell)));
        end
    end
end

end

%%
%Extract the indexes of the entries in the BRENDA data that meet the conditions specified by the search criteria
function EC_indexes = extract_indexes(EC_indCellStringMatches,subs_cell,orgs_cell,subs,...
    substrate,organism, org_index,...
    phylDist)

EC_indexes = EC_indCellStringMatches;%reuse so the string comparisons are not run many times

%If substrate=true then it will extract only the substrates appereances
%indexes in the EC subset from the BRENDA cell array
if substrate
    if (~isempty(EC_indexes)) %optimization
        Subs_indexes = [];
        for l = 1:length(subs)
            if (isempty(subs{l}))
                break;
            end
            Subs_indexes = horzcat(Subs_indexes,EC_indexes(strcmpi(subs(l),...
                subs_cell(EC_indexes))));
        end
        EC_indexes = Subs_indexes;
    end
end

EC_orgs = orgs_cell(EC_indexes);
%If specific organism values are requested looks for all the organism
%repetitions on the subset BRENDA cell array(EC_indexes)
if string(organism) ~= ''
    EC_indexes = EC_indexes(strcmpi(string(organism),EC_orgs));

    %If KEGG code was assigned to the organism (model) then it will look for
    %the Kcat value for the closest organism
elseif org_index~='*' %&& org_index~=''
    KEGG_indexes = [];temp = [];

    %For relating a phyl dist between the modelled organism and the organisms
    %on the BRENDA cell array it should search for a KEGG code for each of
    %these
    for j=1:length(EC_indexes)
        %Assigns a KEGG index for those found on the KEGG struct
        orgs_index = find(strcmpi(orgs_cell(EC_indexes(j)),phylDist.names),1);
        if ~isempty(orgs_index)
            KEGG_indexes = [KEGG_indexes; orgs_index];
            temp         = [temp;EC_indexes(j)];
            %For values related to organisms without KEGG code, then it will
            %look for KEGG code for the first organism with the same genus
        else
            org = orgs_cell{EC_indexes(j)};
            orgGenus = lower(regexprep(org,'\s.*',''));
            if isKey(phylDist.genusHashMap,orgGenus) %annoyingly, this seems to be needed
                matchInd = cell2mat(values(phylDist.genusHashMap,{orgGenus}));
                matches = phylDist.uniqueGenusIndices{matchInd};
                k = matches(1);
                KEGG_indexes = [KEGG_indexes;k];
                temp         = [temp;EC_indexes(j)];
            end
        end
    end
    %Update the EC_indexes cell array
    EC_indexes = temp;
    %Looks for the taxonomically closest organism and saves the index of
    %its appearences in the BRENDA cell
    if ~isempty(EC_indexes)
        distances = phylDist.distMat(org_index,KEGG_indexes);
        EC_indexes = EC_indexes(distances == min(distances));
    end
end
end