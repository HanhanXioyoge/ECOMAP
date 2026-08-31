function prepared = okoPrepareModel(model, biomassRxn, targetRxn, intervals, options)
% okoPrepareModel  Validate an ecModel and construct OKO enzyme-reaction pairs.
% This public helper intentionally contains both legacy and ECOMAP adapters,
% so mapping and unit conversion can be tested without invoking Gurobi.

    if nargin < 4, intervals = []; end
    if nargin < 5, options = struct(); end
    required = {'S','rxns','mets','lb','ub','c'};
    for i = 1:numel(required)
        if ~isfield(model, required{i})
            error('okoPrepareModel:MissingField', 'Model is missing field "%s".', required{i});
        end
    end
    model.b = fieldOr(model, 'b', zeros(size(model.S,1),1));
    biomassIdx = resolveReaction(model, biomassRxn, 'biomass');
    targetIdx = resolveReaction(model, targetRxn, 'target');

    profile = lower(char(fieldOr(options, 'Profile', 'auto')));
    if strcmp(profile, 'auto')
        if isfield(model, 'enzymeConstraints') && ...
                isfield(model.enzymeConstraints, 'ecModeltype') && ...
                strcmpi(model.enzymeConstraints.ecModeltype, 'integrated')
            profile = 'integrated';
        else
            profile = 'legacy';
        end
    end
    if contains(profile, 'integrated')
        [pairs, complexGroups, enzymeVars] = integratedPairs(model);
        enzymeVarSign = -ones(numel(enzymeVars),1);
        canonicalProfile = 'integrated';
    elseif contains(profile, 'legacy')
        [pairs, complexGroups, enzymeVars] = legacyPairs(model);
        enzymeVarSign = ones(numel(enzymeVars),1);
        canonicalProfile = profile;
    else
        error('okoPrepareModel:UnknownProfile', 'Unknown OKO profile "%s".', profile);
    end
    if isempty(pairs)
        error('okoPrepareModel:NoKcatPairs', 'No enzyme-reaction kcat pairs were found.');
    end

    [pairs, intervalInfo] = attachIntervals(pairs, intervals, model, canonicalProfile);
    prepared = struct('model', model, 'profile', canonicalProfile, ...
        'biomassRxn', char(biomassRxn), 'targetRxn', char(targetRxn), ...
        'biomassIdx', biomassIdx, 'targetIdx', targetIdx, ...
        'pairs', pairs, 'complexGroups', {complexGroups}, ...
        'enzymeVarIdx', enzymeVars(:), 'enzymeVarSign', enzymeVarSign, ...
        'intervalInfo', intervalInfo);
end

function [pairs, groups, enzymeVars] = legacyPairs(model)
    if ~isfield(model, 'metNames') || ~isfield(model, 'enzymes')
        error('okoPrepareModel:LegacyFields', ...
            'Legacy profile requires model.metNames and model.enzymes.');
    end
    names = cellstr(string(model.metNames));
    prot = find(contains(names, 'prot_'));
    pmet = find(contains(names, 'pmet_'));
    if isempty(prot) && isempty(pmet)
        error('okoPrepareModel:LegacyLayout', 'Legacy protein metabolite rows were not found.');
    end
    firstEnzMet = min([prot(:); pmet(:)]);
    firstUsage = size(model.S,2) - numel(model.enzymes);
    metabolicLast = firstUsage - 1;
    proteinRows = prot(prot < size(model.S,1)); % exclude prot_pool if it is last
    [rr, cc] = find(model.S(proteinRows, 1:metabolicLast) < 0);
    metIdx = proteinRows(rr);
    coeff = full(model.S(sub2ind(size(model.S), metIdx, cc)));
    keep = metIdx >= firstEnzMet & coeff < 0 & isfinite(coeff);
    metIdx = metIdx(keep); cc = cc(keep); coeff = coeff(keep);
    enzymeID = erase(string(model.metNames(metIdx)), ["prot_","pmet_"]);
    % Legacy coefficients are -1/(3600*kcat); keep public kcat units in s^-1.
    pairs = makePairTable(model, metIdx, cc, coeff, ones(size(coeff))/3600, enzymeID);
    groups = groupsByReaction(cc);
    enzymeVars = (firstUsage:firstUsage + numel(model.enzymes) - 1)';
    enzymeVars = enzymeVars(enzymeVars >= 1 & enzymeVars <= size(model.S,2));
end

function [pairs, groups, enzymeVars] = integratedPairs(model)
    if ~isfield(model, 'enzymeConstraints')
        error('okoPrepareModel:IntegratedFields', 'Integrated profile requires enzymeConstraints.');
    end
    ec = model.enzymeConstraints;
    req = {'rxns','enzymes','rxnEnzMat','kcat'};
    for i = 1:numel(req)
        if ~isfield(ec, req{i})
            error('okoPrepareModel:IntegratedFields', ...
                'enzymeConstraints is missing field "%s".', req{i});
        end
    end
    if isfield(ec, 'ecModeltype') && ~strcmpi(ec.ecModeltype, 'integrated')
        error('okoPrepareModel:UnsupportedModelType', ...
            'Exact OKO requires an integrated ecModel; received "%s".', ec.ecModeltype);
    end
    [ecRow, enzCol] = find(ec.rxnEnzMat ~= 0);
    [foundRxn, rxnIdx] = ismember(ec.rxns(ecRow), model.rxns);
    [foundMet, metIdx] = ismember(strcat('prot_', ec.enzymes(enzCol)), model.mets);
    valid = foundRxn & foundMet & ec.kcat(ecRow) > 0;
    ecRow = ecRow(valid); enzCol = enzCol(valid);
    rxnIdx = rxnIdx(valid); metIdx = metIdx(valid);
    coeff = full(model.S(sub2ind(size(model.S), metIdx, rxnIdx)));
    valid = coeff < 0 & isfinite(coeff);
    ecRow = ecRow(valid); enzCol = enzCol(valid);
    rxnIdx = rxnIdx(valid); metIdx = metIdx(valid); coeff = coeff(valid);
    kcat = double(ec.kcat(ecRow));
    factor = -coeff .* kcat; % MW*subunits/3600 for ECOMAP integrated models
    pairs = makePairTable(model, metIdx, rxnIdx, coeff, factor, string(ec.enzymes(enzCol)));
    pairs.ecRow = ecRow;
    pairs.enzymeColumn = enzCol;
    groups = groupsByReaction(rxnIdx);
    usageIDs = strcat('usage_prot_', ec.enzymes(:));
    [tf, enzymeVars] = ismember(usageIDs, model.rxns);
    if ~any(tf)
        usageIDs = strcat('draw_prot_', ec.enzymes(:));
        [tf, enzymeVars] = ismember(usageIDs, model.rxns);
    end
    enzymeVars = enzymeVars(tf);
    if isempty(enzymeVars)
        error('okoPrepareModel:NoEnzymeVariables', ...
            'Integrated enzyme usage reactions (usage_prot_* or draw_prot_*) were not found.');
    end
end

function pairs = makePairTable(model, metIdx, rxnIdx, coeff, factor, enzymeID)
    kcat = factor ./ (-coeff);
    rxnName = string(model.rxns(rxnIdx));
    if isfield(model, 'rxnNames')
        descriptive = string(model.rxnNames(rxnIdx));
    else
        descriptive = rxnName;
    end
    pairs = table((1:numel(rxnIdx))', metIdx(:), rxnIdx(:), enzymeID(:), ...
        rxnName(:), descriptive(:), coeff(:), factor(:), kcat(:), ...
        nan(numel(rxnIdx),1), nan(numel(rxnIdx),1), false(numel(rxnIdx),1), ...
        'VariableNames', {'pairID','metIdx','rxnIdx','enzymeID','rxnID', ...
        'rxnName','coefficient','factor','kcat','kcatMin','kcatMax','hasInterval'});
end

function groups = groupsByReaction(rxnIdx)
    groups = {};
    ids = unique(rxnIdx(:)');
    for i = 1:numel(ids)
        members = find(rxnIdx == ids(i));
        if numel(members) > 1, groups{end+1,1} = members(:)'; end %#ok<AGROW>
    end
end

function [pairs, info] = attachIntervals(pairs, intervals, model, profile)
    info = struct('provided', false, 'rows', 0, 'matchedPairs', 0, 'unmatchedPairs', height(pairs));
    if isempty(intervals), return; end
    if ischar(intervals) || (isstring(intervals) && isscalar(intervals))
        intervals = readtable(char(intervals), 'VariableNamingRule', 'preserve');
    end
    if ~istable(intervals), error('okoPrepareModel:IntervalType', 'Intervals must be a table or CSV path.'); end
    names = lower(regexprep(string(intervals.Properties.VariableNames), '[^a-z0-9]', ''));
    rxnCol = find(ismember(names, ["rxn","reaction","reactionid","rxnid","rxnname"]),1);
    enzCol = find(ismember(names, ["uniprot","enzyme","enzymeid","protein","proteinid"]),1);
    minCol = find(ismember(names, ["min","minimum","kcatmin","minkcat"]),1);
    maxCol = find(ismember(names, ["max","maximum","kcatmax","maxkcat"]),1);
    if isempty(rxnCol) || isempty(enzCol) || isempty(minCol) || isempty(maxCol)
        error('okoPrepareModel:IntervalColumns', ...
            'Intervals require reaction, enzyme/UniProt, min and max columns.');
    end
    irxn = string(intervals{:,rxnCol}); ienz = string(intervals{:,enzCol});
    imin = numericColumn(intervals{:,minCol}); imax = numericColumn(intervals{:,maxCol});
    if contains(profile, 'legacy')
        [irxn, ienz, imin, imax] = expandLegacyArmRows(irxn, ienz, imin, imax, model);
    end
    for p = 1:height(pairs)
        rxnMatch = strcmpi(strtrim(irxn), pairs.rxnID(p)) | strcmpi(strtrim(irxn), pairs.rxnName(p));
        enzMatch = contains(upper(ienz), upper(pairs.enzymeID(p)));
        hit = find(rxnMatch & enzMatch & isfinite(imin) & isfinite(imax) & imin > 0 & imax > 0);
        if ~isempty(hit)
            pairs.kcatMin(p) = min([pairs.kcat(p); imin(hit); imax(hit)]);
            pairs.kcatMax(p) = max([pairs.kcat(p); imin(hit); imax(hit)]);
            pairs.hasInterval(p) = true;
        end
    end
    info.provided = true; info.rows = height(intervals);
    info.matchedPairs = nnz(pairs.hasInterval); info.unmatchedPairs = nnz(~pairs.hasInterval);
end

function [rxns, enzymes, mins, maxs] = expandLegacyArmRows(rxns, enzymes, mins, maxs, model)
    % Released DLKcat CSVs use reaction names for E. coli and arm reaction
    % IDs for both organisms. Resolve names first, then expand arm entries
    % to the enzyme-specific No* reaction exactly as the legacy scripts do.
    if isfield(model, 'rxnNames')
        for i = 1:numel(rxns)
            hit = find(strcmp(model.rxnNames, rxns(i)), 1);
            if ~isempty(hit), rxns(i) = string(model.rxns{hit}); end
        end
    end
    outR = strings(0,1); outE = strings(0,1); outMin = zeros(0,1); outMax = zeros(0,1);
    for i = 1:numel(rxns)
        rxn = char(rxns(i));
        if startsWith(rxn, 'arm_') && isfield(model, 'grRules')
            root = [rxn(5:end) 'No'];
            candidates = find(startsWith(model.rxns, root) & ~strcmp(model.rxns, rxn));
            gene = '';
            if isfield(model, 'enzymes') && isfield(model, 'enzGenes')
                enzHit = find(contains(model.enzymes, char(enzymes(i))), 1);
                if ~isempty(enzHit), gene = model.enzGenes{enzHit}; end
            end
            if ~isempty(gene)
                candidates = candidates(contains(model.grRules(candidates), gene));
            end
            for j = 1:numel(candidates)
                outR(end+1,1)=string(model.rxns{candidates(j)}); outE(end+1,1)=enzymes(i); %#ok<AGROW>
                outMin(end+1,1)=mins(i); outMax(end+1,1)=maxs(i); %#ok<AGROW>
            end
        else
            outR(end+1,1)=rxns(i); outE(end+1,1)=enzymes(i); %#ok<AGROW>
            outMin(end+1,1)=mins(i); outMax(end+1,1)=maxs(i); %#ok<AGROW>
        end
    end
    rxns=outR; enzymes=outE; mins=outMin; maxs=outMax;
end

function x = numericColumn(x)
    if isnumeric(x), x = double(x); else, x = str2double(string(x)); end
    x = x(:);
end

function idx = resolveReaction(model, id, role)
    idx = find(strcmp(model.rxns, char(id)), 1);
    if isempty(idx), error('okoPrepareModel:ReactionNotFound', '%s reaction "%s" was not found.', role, char(id)); end
end

function value = fieldOr(s, name, fallback)
    if isfield(s, name) && ~isempty(s.(name)), value = s.(name); else, value = fallback; end
end
