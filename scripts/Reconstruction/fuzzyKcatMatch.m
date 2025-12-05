function kcatList = fuzzyKcatMatch(model, ecRxns, parameters, forceWClvl)
% fuzzyKcatMatch (ECOMAP)
%   Adapted from GECKO's fuzzyKcatMatching. This function searches BRENDA
%   for kcat values given EC numbers and substrates from an enzyme-constrained
%   model and returns a per-reaction kcat suggestion. If no exact match is found,
%   it progressively relaxes specificity by: (i) switching organism within
%   phylogenetically close species; (ii) allowing different substrates; (iii)
%   computing kcat from specific activity (SA*MW); and (iv) introducing wildcards
%   in the EC number (1.1.1.- / 1.1.-.- / 1.-.-.-).
%
%   Origin of this implementation:
%     - This ECOMAP function is adapted from GECKO's fuzzyKcatMatching (GECKO 3).
%       Minor structural changes were made to fit ECOMAP's model fields
%       (model.enzymeConstraints.*) and parameter access patterns.
%
% Inputs
%   model        : ECOMAP ecModel-like structure. Requires fields:
%                   - model.rxns, model.metNames, model.S
%                   - model.enzymeConstraints.rxns, .eccodes, .ecModeltype
%   ecRxns       : logical vector or numeric indices selecting the subset of
%                   model.enzymeConstraints.rxns to process. If empty, use all.
%   parameters   : struct with fields
%                   - parameters.org_name : model organism name (string)
%                   - (paths) used by loadBRENDA()/KEGG_struct via basePath
%   forceWClvl   : integer >= 0. Force a minimum wildcard level for EC numbers.
%                   0 = none (exact EC); 1..3 progressively replace last EC levels
%                   by '-' before searching.
%
% Output
%   kcatList     : struct array (one row per processed ec reaction):
%                   - source      : 'brenda'
%                   - rxns        : reaction IDs (from model.enzymeConstraints.rxns)
%                   - substrates  : cell-of-cellstr; model substrates per reaction
%                   - kcats       : proposed kcat (/s). 0 if not found
%                   - eccodes     : EC string(s) used to query BRENDA (with wildcards if applied)
%                   - wildcardLvl : number of EC wildcards used (0..3, NaN if none)
%                   - origin      : specificity tier (1..6, NaN if none)
%                       1: same organism, same substrate, kcat
%                       2: any organism, same substrate, kcat
%                       3: same organism, any substrate, kcat
%                       4: same organism, SA*MW
%                       5: any organism, any substrate, kcat
%                       6: any organism, SA*MW
%
% Notes
%   - When wildcards are used, substrate specificity (tiers 1–2) is ignored by design.
%   - SA*MW uses BRENDA “max_SA” and “max_MW” tables (see loadBRENDA).
%
% Usage
%   kcatList = fuzzyKcatMatch(model, ecRxns, parameters, forceWClvl)

    % ---------------------- Inputs & defaults ----------------------
    if nargin < 2 || isempty(ecRxns)
        ecRxns = true(numel(model.enzymeConstraints.rxns), 1);
    elseif isnumeric(ecRxns)
        ecRxnsVec = false(numel(model.enzymeConstraints.rxns), 1);
        ecRxnsVec(ecRxns) = true;
        ecRxns = ecRxnsVec;
    end
    ecRxns = find(ecRxns); % convert to indices

    if nargin < 3 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end

    if nargin < 4 || isempty(forceWClvl)
        forceWClvl = 0;
    end

    if ~isfield(model.enzymeConstraints, 'eccodes')
        error('No EC codes found in model.enzymeConstraints.eccodes. Run getECnumber() first.');
    end

    % ---------------------- Collect ECs & substrates ----------------------
    eccodes      = model.enzymeConstraints.eccodes(ecRxns);
    substrates   = cell(numel(ecRxns), 1);
    substrCoeffs = cell(numel(ecRxns), 1);

    % Map enzymeConstraints rxn IDs to base model.rxns (handle ECOMAP basic prefix)
    switch model.enzymeConstraints.ecModeltype
        case {'isozyme','integrated'}
            rxnNames = model.enzymeConstraints.rxns;
        case 'basic'
            rxnNames = extractAfter(model.enzymeConstraints.rxns, 4); % strip 3-digit '_' prefixes e.g., '001_'
        otherwise
            rxnNames = model.enzymeConstraints.rxns;
    end

    [~, originalRxns] = ismember(rxnNames(ecRxns), model.rxns);
    for i = 1:length(ecRxns)
        sel = find(model.S(:, originalRxns(i)) < 0);      % substrates: negative stoichiometry
        substrates{i}   = model.metNames(sel);
        substrCoeffs{i} = -model.S(sel, originalRxns(i)); % positive coefficients
    end

    % ---------------------- Load BRENDA & KEGG resources ------------------
    basePath = fullfile(findECOMAProot, 'scripts', 'database');
    [KCATcell, SAcell] = loadBRENDA(basePath);

    % KEGG phylogenetic distances (for nearest-organism fallback)
    PhylDistStructPath = fullfile(basePath, 'keggPhylDist.mat');
    phylDistStruct = KEGG_struct(PhylDistStructPath);

    % Map model organism name to KEGG index (with genus fallback)
    org_name  = parameters.org_name;
    org_index = find_inKEGG(org_name, phylDistStruct.names);

    % Build genus indices/hashes for genus-level fallback when KEGG code is missing
    phylDistStruct.genus = lower(regexprep(phylDistStruct.names, '\s.*', ''));
    phylDistStruct.uniqueGenusList   = unique(phylDistStruct.genus);
    phylDistStruct.genusHashMap      = containers.Map(phylDistStruct.uniqueGenusList, 1:length(phylDistStruct.uniqueGenusList));
    phylDistStruct.uniqueGenusIndices = cell(length(phylDistStruct.uniqueGenusList), 1);
    for i = 1:length(phylDistStruct.genus)
        matchInd = cell2mat(values(phylDistStruct.genusHashMap, phylDistStruct.genus(i)));
        phylDistStruct.uniqueGenusIndices{matchInd} = [phylDistStruct.uniqueGenusIndices{matchInd}; i];
    end

    % ---------------------- Output buffers & stats ------------------------
    kcats = zeros(length(eccodes), 1);
    mM    = length(eccodes);

    % internal info (used to build wildcardLvl and origin; not returned)
    kcatInfo.info.org_s   = zeros(mM,1);
    kcatInfo.info.rest_s  = zeros(mM,1);
    kcatInfo.info.org_ns  = zeros(mM,1);
    kcatInfo.info.rest_ns = zeros(mM,1);
    kcatInfo.info.org_sa  = zeros(mM,1);
    kcatInfo.info.rest_sa = zeros(mM,1);
    kcatInfo.info.wcLevel = NaN(mM,1);
    kcatInfo.stats.queries = 0;
    kcatInfo.stats.org_s   = 0;
    kcatInfo.stats.rest_s  = 0;
    kcatInfo.stats.org_ns  = 0;
    kcatInfo.stats.rest_ns = 0;
    kcatInfo.stats.org_sa  = 0;
    kcatInfo.stats.rest_sa = 0;
    kcatInfo.stats.wc0 = 0; kcatInfo.stats.wc1 = 0; kcatInfo.stats.wc2 = 0;
    kcatInfo.stats.wc3 = 0; kcatInfo.stats.wc4 = 0;
    kcatInfo.stats.matrix = zeros(6,5); % origin x (wildcards 0..4)

    % ---------------------- EC lookup acceleration ------------------------
    % Build an index for unique EC strings to avoid repeated scans
    [ECIndexIds, ~, ic] = unique(KCATcell{1});
    EcIndexIndices = cell(length(ECIndexIds), 1);
    for i = 1:length(EcIndexIndices)
        EcIndexIndices{i} = find(ic == i).';
    end

    % ---------------------- Enforce initial wildcard level ----------------
    % Apply forceWClvl: progressively replace the last EC levels by '-' before searching
    tmpLvl = forceWClvl;
    while tmpLvl > 0
        eccodes = regexprep(eccodes, '(.)*(\.\d+)(\.-)*$', '$1\.-$3');
        tmpLvl = tmpLvl - 1;
    end
    if forceWClvl == 1
        eccodes = regexprep(eccodes, '.*', '-\.-\.-\.-');
    end

    % ---------------------- Main matching loop ----------------------------
    for i = 1:mM
        EC = eccodes{i};
        if ~isempty(EC)
            EC = strsplit(EC, ';'); % support multiple ECs per reaction
            if ~isempty(substrates{i})
                [kcats(i), kcatInfo.info, kcatInfo.stats] = iterativeMatch( ...
                    EC, substrates{i}, substrCoeffs{i}, i, KCATcell, ...
                    kcatInfo.info, kcatInfo.stats, org_name, ...
                    phylDistStruct, org_index, SAcell, ECIndexIds, EcIndexIndices);
            end
        end
    end

    % ---------------------- Assemble kcatList -----------------------------
    kcatList.source      = 'brenda';
    kcatList.rxns        = model.enzymeConstraints.rxns(ecRxns);
    kcatList.substrates  = substrates;
    kcatList.kcats       = kcats;
    kcatList.eccodes     = eccodes;
    kcatList.wildcardLvl = kcatInfo.info.wcLevel;

    % Derive origin vector (1..6) from internal flags:
    kcatList.origin = NaN(numel(model.enzymeConstraints.rxns(ecRxns)), 1);
    originFlags = [kcatInfo.info.org_s, kcatInfo.info.rest_s, ...
                   kcatInfo.info.org_ns, kcatInfo.info.rest_ns, ...
                   kcatInfo.info.org_sa, kcatInfo.info.rest_sa];
    for i = 1:6
        kcatList.origin(find(originFlags(:, i))) = i; %#ok<FNDSB>
    end
end

% =========================================================================
function [kcat,dir,tot] = iterativeMatch(EC, subs, substrCoeff, i, KCATcell, dir, tot, ...
                                         name, phylDist, org_index, SAcell, ECIndexIds, EcIndexIndices)
% iterativeMatch
%   Try a given EC (or list of ECs) with progressively more wildcards until
%   a match is found; if multiple ECs match, choose the one with:
%     (1) minimal wildcard level -> (2) best origin tier -> (3) maximal kcat.

    kcat    = zeros(size(EC));
    origin  = zeros(size(EC));
    matches = zeros(size(EC));
    wc_num  = ones(size(EC)) * 1000; % large sentinel

    for k = 1:length(EC)
        success = false;
        while ~success
            [kcat(k), origin(k), matches(k)] = mainMatch(EC{k}, subs, substrCoeff, KCATcell, ...
                name, phylDist, org_index, SAcell, ECIndexIds, EcIndexIndices);

            if origin(k) > 0
                success   = true;
                wc_num(k) = sum(EC{k} == '-');
            else
                % escalate wildcard: a.b.c.d -> a.b.c.- -> a.b.-.- -> a.-.-.-
                dot_pos  = [2 strfind(EC{k}, '.')];
                wild_num = sum(EC{k} == '-');
                wc_text  = '-.-.-.-';
                EC{k}    = [EC{k}(1:dot_pos(4 - wild_num)) wc_text(1:2*wild_num + 1)];
            end
        end
    end

    if sum(origin) > 0
        % Choose among multiple ECs: least wildcards -> best origin -> max kcat
        best_pos   = (wc_num == min(wc_num));
        new_origin = origin(best_pos);
        best_pos   = (origin == min(new_origin(new_origin ~= 0)));
        max_pos    = find(kcat == max(kcat(best_pos)));
        wc_num     = wc_num(max_pos(1));
        origin     = origin(max_pos(1));
        matches    = matches(max_pos(1));
        kcat       = kcat(max_pos(1));

        % Record per-reaction info & global stats
        dir.org_s(i)   = matches * (origin == 1);
        dir.rest_s(i)  = matches * (origin == 2);
        dir.org_ns(i)  = matches * (origin == 3);
        dir.org_sa(i)  = matches * (origin == 4);
        dir.rest_ns(i) = matches * (origin == 5);
        dir.rest_sa(i) = matches * (origin == 6);
        dir.wcLevel(i) = wc_num;

        tot.org_s   = tot.org_s   + (origin == 1);
        tot.rest_s  = tot.rest_s  + (origin == 2);
        tot.org_ns  = tot.org_ns  + (origin == 3);
        tot.org_sa  = tot.org_sa  + (origin == 4);
        tot.rest_ns = tot.rest_ns + (origin == 5);
        tot.rest_sa = tot.rest_sa + (origin == 6);

        tot.wc0 = tot.wc0 + (wc_num == 0);
        tot.wc1 = tot.wc1 + (wc_num == 1);
        tot.wc2 = tot.wc2 + (wc_num == 2);
        tot.wc3 = tot.wc3 + (wc_num == 3);
        tot.wc4 = tot.wc4 + (wc_num == 4);

        tot.queries = tot.queries + 1;
        tot.matrix(origin, wc_num + 1) = tot.matrix(origin, wc_num + 1) + 1;
    end
end

% =========================================================================
function [kcat, origin, matches] = mainMatch(EC, subs, substrCoeff, KCATcell, ...
                                             name, phylDist, org_index, SAcell, ECIndexIds, EcIndexIndices)
% mainMatch
%   For a single EC string (possibly with wildcards), try the matching tiers
%   in priority order:
%     1) same organism + same substrate, kcat
%     2) any organism  + same substrate, kcat
%     3) same organism + any substrate, kcat
%     4) same organism + SA*MW
%     5) any organism  + any substrate, kcat
%     6) any organism  + SA*MW

    % Precompute EC string matches once (heavy step)
    wild     = false;
    wild_pos = strfind(EC, '-');
    if ~isempty(wild_pos)
        EC   = EC(1:wild_pos(1) - 1);
        wild = true;
    end
    stringMatchesEC_cell = extract_string_matches(EC, KCATcell{1}, wild, ECIndexIds, EcIndexIndices);

    origin = 0;

    % (1) same organism + same substrate
    [kcat, matches] = matchKcat(EC, subs, substrCoeff, KCATcell, name, true, false, ...
                                phylDist, org_index, SAcell, stringMatchesEC_cell, []);
    if matches > 0 && ~wild
        origin = 1;
    else
        % (2) any organism + same substrate (nearest by KEGG distance)
        [kcat, matches] = matchKcat(EC, subs, substrCoeff, KCATcell, '', true, false, ...
                                    phylDist, org_index, SAcell, stringMatchesEC_cell, []);
        if matches > 0 && ~wild
            origin = 2;
        else
            % (3) same organism + any substrate
            [kcat, matches] = matchKcat(EC, subs, substrCoeff, KCATcell, name, false, false, ...
                                        phylDist, org_index, SAcell, stringMatchesEC_cell, []);
            if matches > 0
                origin = 3;
            else
                % (4) same organism + SA*MW
                stringMatchesSA = extract_string_matches(EC, SAcell{1}, wild, [], []);
                [kcat, matches] = matchKcat(EC, subs, substrCoeff, KCATcell, name, false, true, ...
                                            phylDist, org_index, SAcell, stringMatchesEC_cell, stringMatchesSA);
                if matches > 0
                    origin = 4;
                else
                    % (5) any organism + any substrate
                    [kcat, matches] = matchKcat(EC, subs, substrCoeff, KCATcell, '', false, false, ...
                                                phylDist, org_index, SAcell, stringMatchesEC_cell, stringMatchesSA);
                    if matches > 0
                        origin = 5;
                    else
                        % (6) any organism + SA*MW
                        [kcat, matches] = matchKcat(EC, subs, substrCoeff, KCATcell, '', false, true, ...
                                                    phylDist, org_index, SAcell, stringMatchesEC_cell, stringMatchesSA);
                        if matches > 0
                            origin = 6;
                        end
                    end
                end
            end
        end
    end
end

% =========================================================================
function [kcat, matches] = matchKcat(EC, subs, substrCoeff, KCATcell, organism, ...
                                     substrate, SA, phylDist, ...
                                     org_index, SAcell, KCATcellMatches, SAcellMatches)
% matchKcat
%   Select BRENDA entries under the current constraints and return a single
%   kcat value. If multiple values exist, choose the maximum. Upper-bound
%   kcat at 1e7 (/s) as a diffusion-limit safeguard (Bar-Even et al. 2011).

    kcat    = [];
    matches = 0;

    if SA
        % SA path: use SAcell matches (EC+organism via phylDist if needed)
        EC_indexes = extract_indexes(SAcellMatches, [], SAcell{2}, subs, substrate, ...
                                     organism, org_index, phylDist);
        kcat      = SAcell{3}(EC_indexes);   %#ok<NASGU> % (already SA*MW in loadBRENDA)
        % org_cell = SAcell{2}(EC_indexes);  %#ok<NASGU>
        % MW_BRENDA= SAcell{4}(EC_indexes);  %#ok<NASGU>
    else
        % direct kcat path: use KCATcell matches (EC+substrate+organism)
        EC_indexes = extract_indexes(KCATcellMatches, KCATcell{2}, KCATcell{3}, ...
                                     subs, substrate, organism, org_index, phylDist);
        if substrate
            % enforce substrate-name equality; normalize by minimal stoich coeff
            for j = 1:length(EC_indexes)
                idx = EC_indexes(j);
                for k = 1:length(subs)
                    if isempty(subs{k}), break; end
                    if ~isempty(subs{k}) && strcmpi(subs{k}, KCATcell{2}(idx))
                        if KCATcell{4}(idx) > 0
                            coeff  = min(substrCoeff);
                            kCatTmp = KCATcell{4}(idx);
                            kcat    = [kcat; kCatTmp / coeff]; %#ok<AGROW>
                        end
                    end
                end
            end
        else
            kcat = KCATcell{4}(EC_indexes);
        end
    end

    % Aggregate to a single value (max) and cap at 1e7
    if isempty(kcat)
        kcat = 0;
    else
        matches = length(kcat);
        kcat    = max(kcat);
    end
    if kcat > 1e7
        kcat = 1e7;
    end
end

% =========================================================================
function EC_indexes = extract_string_matches(EC, EC_cell, wild, ECIndexIds, EcIndexIndices)
% extract_string_matches
%   Return row indices for EC string matches. If wild==true, match by prefix
%   (startsWith) to emulate hierarchical EC wildcards; otherwise, exact match.
    EC_indexes = [];
    if wild
        if ~isempty(ECIndexIds)
            X = find(contains(ECIndexIds, EC));
            for j = 1:length(X)
                EC_indexes = [EC_indexes, EcIndexIndices{X(j)}]; %#ok<AGROW>
            end
        else
            for j = 1:length(EC_cell)
                if strfind(EC_cell{j}, EC) == 1 %#ok<STRIFCND>
                    EC_indexes = [EC_indexes, j]; %#ok<AGROW>
                end
            end
        end
    else
        if ~isempty(ECIndexIds)
            mtch = find(strcmpi(EC, ECIndexIds));
            if ~isempty(mtch)
                EC_indexes = EcIndexIndices{mtch};
            end
        else
            if ~isempty(EC_cell)
                EC_indexes = transpose(find(strcmpi(EC, EC_cell)));
            end
        end
    end
end

% =========================================================================
function EC_indexes = extract_indexes(EC_indCellStringMatches, subs_cell, orgs_cell, subs, ...
                                      substrate, organism, org_index, phylDist)
% extract_indexes
%   Filter EC candidates by substrate (exact name) and organism.
%   If organism is empty, select closest organism via KEGG distance matrix.
    EC_indexes = EC_indCellStringMatches;

    % substrate filtering
    if substrate
        if ~isempty(EC_indexes)
            Subs_indexes = [];
            for l = 1:length(subs)
                if isempty(subs{l}), break; end
                Subs_indexes = horzcat(Subs_indexes, EC_indexes(strcmpi(subs(l), subs_cell(EC_indexes)))); %#ok<AGROW>
            end
            EC_indexes = Subs_indexes;
        end
    end

    if isempty(orgs_cell) || isempty(EC_indexes)
        return;
    end

    EC_orgs = orgs_cell(EC_indexes);

    % same-organism filter
    if string(organism) ~= ""
        EC_indexes = EC_indexes(strcmpi(string(organism), EC_orgs));

    % nearest-organism by KEGG distance (if model organism mapped)
    elseif org_index ~= '*'
        KEGG_indexes = []; temp = [];
        for j = 1:length(EC_indexes)
            orgs_index = find(strcmpi(orgs_cell(EC_indexes(j)), phylDist.names), 1);
            if ~isempty(orgs_index)
                KEGG_indexes = [KEGG_indexes; orgs_index]; %#ok<AGROW>
                temp         = [temp; EC_indexes(j)];      %#ok<AGROW>
            else
                % genus fallback when KEGG code not available
                org      = orgs_cell{EC_indexes(j)};
                orgGenus = lower(regexprep(org, '\s.*', ''));
                if isKey(phylDist.genusHashMap, orgGenus)
                    matchInd = cell2mat(values(phylDist.genusHashMap, {orgGenus}));
                    matches  = phylDist.uniqueGenusIndices{matchInd};
                    k        = matches(1);
                    KEGG_indexes = [KEGG_indexes; k];      %#ok<AGROW>
                    temp         = [temp; EC_indexes(j)]; %#ok<AGROW>
                end
            end
        end
        EC_indexes = temp;
        if ~isempty(EC_indexes)
            distances = phylDist.distMat(org_index, KEGG_indexes);
            EC_indexes = EC_indexes(distances == min(distances));
        end
    end
end

% =========================================================================
function org_index = find_inKEGG(org_name, names)
% find_inKEGG
%   Map organism name to KEGG index; if multiple hits, take first; if none,
%   try genus-only; if still none, return '*'.
    org_index = find(strcmpi(org_name, names));
    if numel(org_index) > 1
        org_index = org_index(1);
    elseif isempty(org_index)
        org_g    = regexprep(org_name, '\s.*', '');
        org_index = find(strcmpi(org_g, names));
        if numel(org_index) > 1
            org_index = org_index(1);
        elseif isempty(org_index)
            org_index = '*';
        end
    end
end

% =========================================================================
function phylDistStruct = KEGG_struct(phylpath)
% KEGG_struct
%   Load KEGG phylogenetic distance structure saved as MAT file.
    load(phylpath, 'phylDistStruct');
    phylDistStruct.ids   = transpose(phylDistStruct.ids);
    phylDistStruct.names = transpose(phylDistStruct.names);
    phylDistStruct.names = regexprep(phylDistStruct.names, '\s*\(.*', '');
end

% =========================================================================
function [KCATcell, SAcell] = loadBRENDA(basePath)
% loadBRENDA
%   Load preprocessed BRENDA “max_KCAT”, “max_SA”, “max_MW” tables from text
%   files located under basePath, and return KCATcell and SAcell structures.
%   Units:
%     - max_SA.txt is scaled 1/60 to convert [umol/min/mg] -> [mmol/s/g]
%     - max_MW.txt is scaled 1/1000 to convert [g/mol] -> [g/mmol]
%     - SAcell{3} will contain SA*MW (→ /s) for convenience.
%
%   Expected formats (tab-delimited):
%     KCATcell: {EC, SubstrateName, Organism, kcat_value, ...}
%     SAcell : {EC, <unused>, Organism, SA_value, MW_value}
%
%   Note: This loader is adapted to the ECOMAP repository layout.

    KCAT_file = fullfile(basePath, 'max_KCAT.txt');
    SA_file   = fullfile(basePath, 'max_SA.txt');
    MW_file   = fullfile(basePath, 'max_MW.txt');

    KCATcell       = openDataFile(KCAT_file, 1);       % kcat already in /s
    scalingFactor  = 1/60;    % [umol/min/mg] -> [mmol/s/g]
    SA              = openDataFile(SA_file, scalingFactor);
    scalingFactor  = 1/1000;  % [g/mol] -> [g/mmol]
    MW              = openDataFile(MW_file, scalingFactor);

    for i = 1:4
        SAcell{i} = [];
    end
    previousEC = []; EC_indexes = [];

    % Build index on MW{1} (EC) to speed SA*MW lookup
    MWECNum          = upper(unique(MW{1}));
    MWECNumIndices   = cell(length(MWECNum), 1);
    MWECNumHashMap   = containers.Map(MWECNum, 1:length(MWECNum));
    for i = 1:length(MW{1})
        matchInd = cell2mat(values(MWECNumHashMap, MW{1}(i)));
        MWECNumIndices{matchInd} = [MWECNumIndices{matchInd}; i];
    end

    % Build SAcell (EC, Organism, SA*MW, MW)
    for i = 1:length(SA{1})
        if ~strcmpi(SA{1}(i), previousEC)
            key = upper(SA{1}(i));
            if isKey(MWECNumHashMap, key)
                matchInd  = cell2mat(values(MWECNumHashMap, key));
                EC_indexes = MWECNumIndices{matchInd};
            else
                EC_indexes = [];
            end
        end
        mwEC{1} = MW{3}(EC_indexes); % organisms
        mwEC{2} = MW{4}(EC_indexes); % MW (g/mmol)

        org_index = find(strcmpi(SA{3}(i), mwEC{1}), 1);
        if ~isempty(org_index)
            SAcell{1} = [SAcell{1}; SA{1}(i)];                    % EC
            SAcell{2} = [SAcell{2}; SA{3}(i)];                    % Organism
            SAcell{3} = [SAcell{3}; SA{4}(i) * mwEC{2}(org_index)]; % SA*MW (/s)
            SAcell{4} = [SAcell{4}; mwEC{2}(org_index)];          % MW (g/mmol)
        end
        previousEC = SA{1}(i);
    end

    % Strip leading 'EC' from EC numbers (e.g., 'EC 1.1.1.1' -> '1.1.1.1')
    if ~isempty(KCATcell{1}), KCATcell{1} = extractAfter(KCATcell{1}, 2); end
    if ~isempty(SAcell{1}),  SAcell{1}  = extractAfter(SAcell{1}, 2);  end

    function data_cell = openDataFile(fileName, scalingFactor)
        fID       = fopen(fileName);
        data_cell = textscan(fID, '%q %q %q %f %q', 'delimiter', '\t');
        fclose(fID);
        data_cell{4} = data_cell{4} * scalingFactor;
        % Clean trailing comments in organism field
        data_cell{3} = regexprep(data_cell{3}, '\/\/.*', '');
    end
end
