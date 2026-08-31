function csvPath = buildCrossSpeciesDlInput(seeds, candidates, outputDir, predictor)
%BUILDCROSSSPECIESDLINPUT Expand seeds × candidates × substrates into DL input CSV.
%
%   For each (seed) in the input:
%     For each homolog in candidates matching this seed:
%       For each substrate in the seed's Substrates:
%         Write one CSV row with:
%           ReactionName  <- from seed (verbatim, NOT cleaned)
%           Organism      <- from homolog
%           GeneID        <- from homolog (may be empty)
%           ProteinID     <- homolog accession used by the predictor
%           EC Number     <- from seed
%           MetaNetXID    <- from seed substrate
%           Substrate     <- from seed substrate
%           SMILES        <- from seed substrate
%           InChIKey      <- from seed substrate
%           sequence      <- from homolog
%
%   Output CSV has the standard 10 columns (DLKcat/UniKP) or 11 columns (CatPred
%   with pdbpath). Compatible with ExecutePrediction.m.
%
%   csvPath = buildCrossSpeciesDlInput(seeds, candidates, outputDir, predictor)
%
%   Inputs:
%     seeds       - table from buildOkoReactionSeeds with columns:
%                   ReactionName, OriginalProteinID, ECNumber,
%                   Substrates (cellstr), SubstrateIdentities (cellstr)
%     candidates  - table from fetchCrossSpeciesEnzymes with columns:
%                   ReactionName, OriginalProteinID, accession, organism,
%                   sequence, ...
%     outputDir   - directory to write the CSV into
%     predictor   - 'DLKcat' | 'UniKP' | 'CatPred'
%
%   Output:
%     csvPath - full path to written CSV
%                (UniKP_OKOplus_input.csv or CatPred_OKOplus_input.csv)

    if ~ismember(predictor, {'DLKcat','UniKP','CatPred'})
        error('buildCrossSpeciesDlInput:invalidPredictor', ...
              'predictor must be DLKcat, UniKP or CatPred');
    end

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Build lookup: seed identity → substrate metadata
    % We need Substrate, MetaNetXID, SMILES, InChIKey from the DL input.
    % These are NOT in seeds (which has only Substrates names + identities).
    % We rebuild them from SubstrateIdentities and Substrates lists.

    out = table();
    out.ReactionName      = cell(0, 1);
    out.Organism          = cell(0, 1);
    out.GeneID            = cell(0, 1);
    out.ProteinID         = cell(0, 1);
    out.ECNumber          = cell(0, 1);
    out.MetaNetXID        = cell(0, 1);
    out.Substrate         = cell(0, 1);
    out.SMILES            = cell(0, 1);
    out.InChIKey          = cell(0, 1);
    out.sequence          = cell(0, 1);

    if strcmp(predictor, 'CatPred')
        out.pdbpath = cell(0, 1);
    end

    for i = 1:height(seeds)
        seed = seeds(i, :);
        rxnId    = seed.ReactionName{1};
        origProt = seed.OriginalProteinID{1};
        ecNum    = seed.ECNumber{1};

        % Candidates matching this seed
        candMask = strcmp(candidates.ReactionName, rxnId);
        if ismember('OriginalProteinID', candidates.Properties.VariableNames)
            candMask = candMask & strcmp(candidates.OriginalProteinID, origProt);
        end
        candRows = candidates(candMask, :);
        if isempty(candRows), continue; end

        substrates     = seed.Substrates{1};         % cellstr of names
        identities     = seed.SubstrateIdentities{1}; % cellstr of preferred identity
        if ismember('MetaNetXIDs', seeds.Properties.VariableNames)
            mnxids = seed.MetaNetXIDs{1};
            smiles = seed.SMILESList{1};
            inchikeys = seed.InChIKeys{1};
        else
            % Backward-compatible fallback for callers constructing legacy
            % seed tables manually.
            mnxids = identities;
            smiles = identities;
            inchikeys = identities;
        end
        nSubstrates    = numel(substrates);

        for j = 1:height(candRows)
            organism = candRows.organism{j};
            geneId   = '';
            if ismember('geneID', candRows.Properties.VariableNames)
                geneId = candRows.geneID{j};
            end
            sequence = candRows.sequence{j};
            homologAccession = candRows.accession{j};

            for s = 1:nSubstrates
                row = struct();
                row.ReactionName = {rxnId};
                row.Organism     = {organism};
                row.GeneID       = {geneId};
                % Predictor rows describe the homolog sequence.  The original
                % ecModel UniProt stays in the seed table and is restored only
                % after the ReactionName-level aggregation.
                row.ProteinID    = {homologAccession};
                row.ECNumber     = {ecNum};
                row.MetaNetXID   = {mnxids{s}};
                row.Substrate    = {substrates{s}};
                row.SMILES       = {smiles{s}};
                row.InChIKey     = {inchikeys{s}};
                row.sequence     = {sequence};

                % Try to recover SMILES / MNXID from candidates table for this organism
                % (if candidates has a SMILES/MNXID column). Otherwise leave as identity.
                % For now, use identity as fallback (works for InChIKey-rich inputs).

                if strcmp(predictor, 'CatPred')
                    % Cross-species homologs do not have a model-provided
                    % structure path. Use CatPred's explicit sentinel value.
                    row.pdbpath = {'none'};
                end

                out = [out; struct2table(row, 'AsArray', true)]; %#ok<AGROW>
            end
        end
    end

    % The same homolog/substrate predictor row can be reached through more
    % than one original ecModel protein. Predict it only once; the complete
    % candidate table restores every (ReactionName, OriginalProteinID) link
    % after prediction. Besides preventing statistical double-counting, this
    % substantially reduces Docker input size and prediction time.
    if ~isempty(out)
        out = unique(out, 'rows', 'stable');
    end

    csvName = [predictor '_OKOplus_input.csv'];
    csvPath = fullfile(outputDir, csvName);
    atomicWriteTable(out, csvPath);
end

function atomicWriteTable(T, targetPath)
    tempPath = [tempname(fileparts(targetPath)) '.csv'];
    cleanup = onCleanup(@()deleteIfPresent(tempPath)); %#ok<NASGU>
    writetable(T, tempPath);
    [ok, message] = movefile(tempPath, targetPath, 'f');
    if ~ok
        error('buildCrossSpeciesDlInput:atomicWriteFailed', '%s', message);
    end
end

function deleteIfPresent(path)
    if exist(path, 'file'), delete(path); end
end
