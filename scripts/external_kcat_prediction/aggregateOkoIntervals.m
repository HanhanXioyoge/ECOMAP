function [aggregated, qc] = aggregateOkoIntervals(predictionsCsv, seeds, predictor, candidates)
%AGGREGATEOKOINTERVALS Build OKO+ intervals per enzyme-reaction pair.
%
%   [aggregated,qc] = aggregateOkoIntervals(csv,seeds,predictor,candidates)
%
%   Predictor ProteinID values are homolog accessions.  The compact mapping
%   (ReactionName, accession, OriginalProteinID) from candidates restores the
%   ecModel enzyme, including legitimate many-to-many mappings. Predictions
%   are de-duplicated before that expansion and finally grouped by
%   (ReactionName, OriginalProteinID), matching the original OKO+ semantics.
%
%   The legacy three-argument form remains available when every reaction in
%   seeds maps to one original protein.

    if nargin < 4, candidates = table(); end
    if ~exist(predictionsCsv, 'file')
        error('aggregateOkoIntervals:missingFile', ...
              'Predictions CSV not found: %s', predictionsCsv);
    end

    kcatCol = predictorKcatColumn(predictor);
    T = readCompactPredictionTable(predictionsCsv, kcatCol);
    if isempty(T)
        [aggregated, qc] = emptyOutputs();
        return;
    end

    requireColumns(T, {'ReactionName','ProteinID',kcatCol}, 'predictions CSV');

    % Keep only compact columns needed for mapping, aggregation and QC. This
    % drops large sequence/SMILES fields before joining a genome-scale file.
    P = table(textVector(T.ReactionName), textVector(T.ProteinID), ...
              toNumericColumn(T.(kcatCol)), optionalText(T, 'Organism'), ...
              substrateIdentity(T), ...
              'VariableNames', {'ReactionName','accession','kcat', ...
                                'Organism','SubstrateKey'});
    P = P(isfinite(P.kcat) & P.kcat > 0 & ...
          P.ReactionName ~= "" & P.accession ~= "", :);
    if isempty(P)
        [aggregated, qc] = emptyOutputs();
        return;
    end

    % Multiple original ecModel proteins can deliberately generate the same
    % predictor row. Predict it once statistically, then fan it back out via
    % the candidate mapping so mean/mode are not biased by duplicate rows.
    P = unique(P, 'rows', 'stable');

    if ~isempty(candidates)
        requireColumns(candidates, ...
            {'ReactionName','OriginalProteinID','accession'}, 'candidate table');
        M = table(textVector(candidates.ReactionName), ...
                  textVector(candidates.accession), ...
                  textVector(candidates.OriginalProteinID), ...
                  'VariableNames', {'ReactionName','accession','OriginalProteinID'});
        M = unique(M(M.ReactionName ~= "" & M.accession ~= "" & ...
                     M.OriginalProteinID ~= "", :), 'rows', 'stable');
        mapped = innerjoin(P, M, 'Keys', {'ReactionName','accession'});
    else
        mapped = legacyMapByReaction(P, seeds);
    end

    if isempty(mapped)
        [aggregated, qc] = emptyOutputs();
        return;
    end

    [G, groupRxn, groupProtein] = findgroups( ...
        mapped.ReactionName, mapped.OriginalProteinID);
    nGroups = numel(groupRxn);
    aggregated = table(cellstr(groupRxn), cellstr(groupProtein), ...
        zeros(nGroups,1), zeros(nGroups,1), zeros(nGroups,1), zeros(nGroups,1), ...
        'VariableNames', {'rxn','uniprot','min','max','mean','mode'});
    qc = table(cellstr(groupRxn), cellstr(groupProtein), ...
        zeros(nGroups,1), zeros(nGroups,1), zeros(nGroups,1), zeros(nGroups,1), ...
        'VariableNames', {'rxn','uniprot','nOrganisms','nProteins', ...
                          'nSubstrates','nPredictions'});

    for g = 1:nGroups
        mask = (G == g);
        k = mapped.kcat(mask);
        aggregated.min(g) = min(k);
        aggregated.max(g) = max(k);
        aggregated.mean(g) = mean(k);
        aggregated.mode(g) = computeModeLog10(k);
        qc.nOrganisms(g) = countNonemptyUnique(mapped.Organism(mask));
        qc.nProteins(g) = countNonemptyUnique(mapped.accession(mask));
        qc.nSubstrates(g) = countNonemptyUnique(mapped.SubstrateKey(mask));
        qc.nPredictions(g) = numel(k);
    end
end

function T = readCompactPredictionTable(csvPath, kcatCol)
% Avoid loading sequence and other large predictor payload columns.
    options = detectImportOptions(csvPath, 'VariableNamingRule', 'preserve');
    required = {'ReactionName','ProteinID',kcatCol};
    missing = setdiff(required, options.VariableNames);
    if ~isempty(missing)
        error('aggregateOkoIntervals:missingColumn', ...
            'predictions CSV missing required columns: %s', ...
            strjoin(missing, ', '));
    end
    wanted = [required, {'Organism','InChIKey','MetaNetXID', ...
                         'Substrate','Substrate_name'}];
    options.SelectedVariableNames = intersect( ...
        wanted, options.VariableNames, 'stable');
    T = readtable(csvPath, options);
end

function mapped = legacyMapByReaction(P, seeds)
% Compatibility path for callers without a candidate mapping.
    requireColumns(seeds, {'ReactionName','OriginalProteinID'}, 'seed table');
    M = unique(table(textVector(seeds.ReactionName), ...
        textVector(seeds.OriginalProteinID), ...
        'VariableNames', {'ReactionName','OriginalProteinID'}), 'rows', 'stable');
    [~, rxn] = findgroups(M.ReactionName);
    counts = splitapply(@(x)numel(unique(x)), M.OriginalProteinID, ...
                        findgroups(M.ReactionName));
    ambiguous = rxn(counts > 1);
    if ~isempty(ambiguous)
        error('aggregateOkoIntervals:missingCandidateMapping', ...
            ['ReactionName "%s" maps to multiple original proteins. Pass ' ...
             'OKOplus_candidates to restore enzyme-reaction pairs.'], ambiguous(1));
    end
    mapped = innerjoin(P, M, 'Keys', 'ReactionName');
end

function col = predictorKcatColumn(predictor)
    switch char(predictor)
        case {'DLKcat','UniKP'}
            col = 'predicted_kcat';
        case 'CatPred'
            col = 'Prediction_(s^(-1))';
        otherwise
            error('aggregateOkoIntervals:invalidPredictor', ...
                'Predictor must be DLKcat, UniKP or CatPred.');
    end
end

function values = optionalText(T, name)
    if ismember(name, T.Properties.VariableNames)
        values = textVector(T.(name));
    else
        values = repmat("", height(T), 1);
    end
end

function key = substrateIdentity(T)
    key = repmat("", height(T), 1);
    aliases = {'InChIKey','MetaNetXID','Substrate','Substrate_name'};
    for i = 1:numel(aliases)
        if ~ismember(aliases{i}, T.Properties.VariableNames), continue; end
        value = textVector(T.(aliases{i}));
        use = key == "" & value ~= "" & ~strcmpi(value, "NA");
        key(use) = value(use);
    end
end

function requireColumns(T, names, source)
    missing = setdiff(names, T.Properties.VariableNames);
    if ~isempty(missing)
        error('aggregateOkoIntervals:missingColumn', ...
            '%s missing required columns: %s', source, strjoin(missing, ', '));
    end
end

function values = textVector(values)
    if ischar(values), values = string(cellstr(values)); else, values = string(values(:)); end
    values(ismissing(values)) = "";
end

function n = countNonemptyUnique(values)
    values = values(values ~= "");
    n = numel(unique(values));
end

function value = computeModeLog10(kcats)
    if numel(kcats) < 2, value = kcats(1); return; end
    logK = log10(kcats);
    try
        [f, xi] = ksdensity(logK);
        [~, idx] = max(f);
        value = 10 ^ xi(idx);
    catch
        value = median(kcats);
    end
end

function v = toNumericColumn(col)
    if isnumeric(col), v = double(col(:)); return; end
    if iscell(col)
        v = nan(numel(col),1);
        for i = 1:numel(col)
            value = col{i};
            if isnumeric(value), v(i) = double(value); else, v(i) = str2double(string(value)); end
        end
        return;
    end
    v = str2double(string(col(:)));
end

function [aggregated, qc] = emptyOutputs()
    aggregated = table(cell(0,1), cell(0,1), zeros(0,1), zeros(0,1), ...
        zeros(0,1), zeros(0,1), ...
        'VariableNames', {'rxn','uniprot','min','max','mean','mode'});
    qc = table(cell(0,1), cell(0,1), zeros(0,1), zeros(0,1), ...
        zeros(0,1), zeros(0,1), ...
        'VariableNames', {'rxn','uniprot','nOrganisms','nProteins', ...
                          'nSubstrates','nPredictions'});
end
