function seeds = buildOkoReactionSeeds(dlInputCsv, predictor)
%BUILDOKOREACTIONSEEDS Extract seed table from existing DL input CSV.
%
%   The seed table encodes the (ReactionName, OriginalProteinID, EC,
%   substrate-set) tuples that OKO+ needs cross-species kcat intervals
%   for. It does NOT enter ExecutePrediction; it is consumed by
%   fetchCrossSpeciesEnzymes and aggregateOkoIntervals.
%
%   seeds = buildOkoReactionSeeds(dlInputCsv)
%   seeds = buildOkoReactionSeeds(dlInputCsv, predictor)
%
%   Inputs:
%     dlInputCsv - path to existing UniKP_input.csv or CatPred_input.csv
%                  (10 columns: ReactionName, Organism, GeneID, ProteinID,
%                   EC Number, MetaNetXID, Substrate, SMILES, InChIKey,
%                   sequence, plus pdbpath for CatPred)
%     predictor  - 'UniKP' | 'CatPred' (default 'UniKP'). Controls which
%                  sequence column is treated as canonical and which
%                  substrate identity column is preferred (InChIKey >
%                  MetaNetXID > SMILES).
%
%   Output:
%     seeds - table with columns:
%       ReactionName       - preserved verbatim from input (NOT cleaned)
%       OriginalProteinID  - target enzyme UniProt from input's ProteinID
%       ECNumber           - EC annotation from input's 'EC Number'
%       Substrates         - cellstr of unique substrate names per rxn
%       SubstrateIdentities- cellstr of preferred identity (InChIKey if
%                            present, else MetaNetXID, else SMILES)
%       nSubstrates        - count of unique substrates per rxn
%
%   See: docs/superpowers/specs/2026-08-09-oko-plus-homolog-interval-pipeline-design.md
%        (refactor: input from existing DL CSV)

    if nargin < 2
        predictor = 'UniKP';
    end

    if ~exist(dlInputCsv, 'file')
        error('buildOkoReactionSeeds:missingFile', ...
              'DL input CSV not found: %s', dlInputCsv);
    end

    T = readtable(dlInputCsv, 'VariableNamingRule', 'preserve');
    if isempty(T)
        seeds = emptySeedsTable();
        return;
    end

    % Required columns per predictor
    required = {'ReactionName', 'ProteinID', 'EC Number', ...
                'Substrate', 'MetaNetXID', 'SMILES', 'InChIKey'};
    missing = setdiff(required, T.Properties.VariableNames);
    if ~isempty(missing)
        error('buildOkoReactionSeeds:missingColumns', ...
              'DL input CSV missing columns: %s', strjoin(missing, ', '));
    end

    % Normalize column names for downstream code: 'EC Number' -> 'ECNumber'
    T.Properties.VariableNames = strrep(T.Properties.VariableNames, ...
                                        'EC Number', 'ECNumber');

    % Group by (ReactionName, ProteinID) — these define an OKO+ target.
    % Use a stable composite key to preserve row order across runs.
    groupKeys = strcat(T.ReactionName, '||', T.ProteinID);
    [uniqueKeys, ~, ic] = unique(groupKeys, 'stable');
    nGroups = numel(uniqueKeys);

    seeds = table( ...
        cell(nGroups, 1), ...  % ReactionName
        cell(nGroups, 1), ...  % OriginalProteinID
        cell(nGroups, 1), ...  % ECNumber
        cell(nGroups, 1), ...  % Substrates (cellstr)
        cell(nGroups, 1), ...  % SubstrateIdentities (cellstr)
        cell(nGroups, 1), ...  % MetaNetXIDs (cellstr)
        cell(nGroups, 1), ...  % SMILESList (cellstr)
        cell(nGroups, 1), ...  % InChIKeys (cellstr)
        zeros(nGroups, 1), ... % nSubstrates
        'VariableNames', {'ReactionName', 'OriginalProteinID', 'ECNumber', ...
                          'Substrates', 'SubstrateIdentities', 'MetaNetXIDs', ...
                          'SMILESList', 'InChIKeys', 'nSubstrates'});

    for g = 1:nGroups
        sepIdx = strfind(uniqueKeys{g}, '||');
        seeds.ReactionName{g}      = uniqueKeys{g}(1:sepIdx-1);
        seeds.OriginalProteinID{g} = uniqueKeys{g}(sepIdx+2:end);
        groupT = T(ic == g, :);

        % EC number: first non-empty value (model uses one per rxn).
        ecVals = groupT.ECNumber;
        seeds.ECNumber{g} = firstNonEmpty(ecVals);

        % Dedupe substrates by preferred identity (InChIKey > MNXID > SMILES).
        % This ensures ATP (same InChIKey) is counted once even if it appears
        % multiple times in the input.
        [uniqueSub, nUnique] = dedupeSubstrates(groupT, predictor);
        seeds.Substrates{g}          = uniqueSub.names;
        seeds.SubstrateIdentities{g} = uniqueSub.identities;
        seeds.MetaNetXIDs{g}         = uniqueSub.mnxids;
        seeds.SMILESList{g}          = uniqueSub.smiles;
        seeds.InChIKeys{g}           = uniqueSub.inchikeys;
        seeds.nSubstrates(g)         = nUnique;
    end
end

function val = firstNonEmpty(vals)
%FIRSTNONEMPTY Return first non-empty entry in a cellstr array.
    val = '';
    if iscell(vals)
        for i = 1:numel(vals)
            v = vals{i};
            if isstring(v), v = char(v); end
            if ~isempty(v)
                val = v;
                return;
            end
        end
    end
end

function [uniqueSub, nUnique] = dedupeSubstrates(groupT, predictor)
%DEDUPESUBSTRATES Dedupe substrates per ReactionName by preferred identity.
%
%   Returns unique substrate names and their preferred-identity strings.
%   Priority: InChIKey > MetaNetXID > SMILES (per user spec).
%   Returns two parallel cellstr arrays and a count.
    if isstring(groupT.Substrate), groupT.Substrate = cellstr(groupT.Substrate); end
    nRows = height(groupT);

    % Build identity preference score (higher = more preferred)
    identity = repmat({''}, nRows, 1);
    for i = 1:nRows
        if isstring(groupT.InChIKey(i)), groupT.InChIKey(i) = cellstr(groupT.InChIKey(i)); end
        if isstring(groupT.MetaNetXID(i)), groupT.MetaNetXID(i) = cellstr(groupT.MetaNetXID(i)); end
        if isstring(groupT.SMILES(i)), groupT.SMILES(i) = cellstr(groupT.SMILES(i)); end

        ik = strtrim(groupT.InChIKey{i});
        if ~isempty(ik) && ~strcmpi(ik, 'NA')
            identity{i} = ik;
            continue;
        end
        mnx = strtrim(groupT.MetaNetXID{i});
        if ~isempty(mnx) && ~strcmpi(mnx, 'NA')
            identity{i} = mnx;
            continue;
        end
        smi = strtrim(groupT.SMILES{i});
        if ~isempty(smi)
            identity{i} = smi;
        end
    end

    % Unique by identity (preserve first-seen order)
    [uniqueIdentities, ia] = unique(identity, 'stable');
    uniqueNames = cellstr(groupT.Substrate(ia));
    uniqueMnx = cellstr(string(groupT.MetaNetXID(ia)));
    uniqueSmiles = cellstr(string(groupT.SMILES(ia)));
    uniqueInchi = cellstr(string(groupT.InChIKey(ia)));

    % Filter out rows with empty identity (defensive: should not happen)
    validMask = ~cellfun(@isempty, uniqueIdentities);
    uniqueNames      = uniqueNames(validMask);
    uniqueIdentities = uniqueIdentities(validMask);
    uniqueMnx        = uniqueMnx(validMask);
    uniqueSmiles     = uniqueSmiles(validMask);
    uniqueInchi      = uniqueInchi(validMask);

    uniqueSub.names      = uniqueNames;
    uniqueSub.identities = uniqueIdentities;
    uniqueSub.mnxids     = uniqueMnx;
    uniqueSub.smiles     = uniqueSmiles;
    uniqueSub.inchikeys  = uniqueInchi;
    nUnique              = numel(uniqueNames);
end

function T = emptySeedsTable()
    T = table(cell(0,1), cell(0,1), cell(0,1), cell(0,1), cell(0,1), ...
              cell(0,1), cell(0,1), cell(0,1), zeros(0,1), ...
              'VariableNames', {'ReactionName', 'OriginalProteinID', 'ECNumber', ...
                                'Substrates', 'SubstrateIdentities', 'MetaNetXIDs', ...
                                'SMILESList', 'InChIKeys', 'nSubstrates'});
end
