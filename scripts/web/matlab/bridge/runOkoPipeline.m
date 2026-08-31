function result = runOkoPipeline(ecModel, predictors, options, parameters)
%RUNOKOPIPELINE  Pure orchestration entry point for the OKO+ homolog pipeline.
%
%   REFACTORED VERSION. Pipeline:
%     1. writeInputFile generates the requested predictor input CSVs.
%     2. buildOkoReactionSeeds extracts one shared seed set.
%     3. fetchCrossSpeciesEnzymes runs ONCE; DLKcat, UniKP and CatPred
%        share the same homolog candidates.
%     4. buildCrossSpeciesDlInput writes UniKP_OKOplus_input.csv /
%        CatPred_OKOplus_input.csv with (rxn × homologs × substrates)
%        expansion.
%     5. ExecutePrediction runs the Docker predictor (or mock).
%     6. aggregateOkoIntervals restores the original enzyme through the
%        candidate mapping and groups by (ReactionName, OriginalProteinID).
%     7. formatOkoPredictionsCsv writes the OKO-compatible 6-column
%        CSV consumed by algOkoPlus.
%
%   result = runOkoPipeline(ecModel, predictors)
%   result = runOkoPipeline(ecModel, predictors, options)
%   result = runOkoPipeline(ecModel, predictors, options, parameters)
%
%   `parameters` is the LAST argument. If empty, the function calls
%   ParameterManager.getParams() to populate it (this is the project's
%   standard convention).
%
%   Inputs:
%     ecModel    - ecModel struct (used to generate the initial DL inputs)
%     params     - struct from ParameterManager.getParams() with at minimum
%                  .reconstructionDir.
%     predictors - cellstr from {'UniKP', 'CatPred'}.
%     options    - optional struct with injection points for testing:
%         .CacheDir    - UniProt cache directory (default <reconstructionDir>/.cache/uniprot)
%         .dockerFn    - @(predictor, outputDir) -> void. Default:
%                        ExecutePrediction({predictor}, outputDir).
%         .Logger      - structured log function. Default: fprintf wrapper.
%         .organismName- char used as CSV filename prefix. Default: ecModel.id.
%         .writeDlInput- @(model, predictor, outputDir) -> csvPath. Default:
%                        writeInputFile(model, predictor, params).
%         .buildSeeds  - @(dlInputCsv, predictor) -> seeds. Default:
%                        buildOkoReactionSeeds.
%         .fetchEnzymes- @(seeds, opts) -> candidates. Default:
%                        fetchCrossSpeciesEnzymes.
%         .buildXlInput- @(seeds, candidates, outputDir, predictor) -> csvPath.
%                        Default: buildCrossSpeciesDlInput.
%         .aggregate   - @(predictionsCsv, seeds, predictor) -> [agg, qc].
%                        Default: aggregateOkoIntervals.
%         .formatCsv   - @(aggregated, csvPath). Default: formatOkoPredictionsCsv.
%
%   Returns struct with:
%     .predictor_csv_paths      - struct mapping predictor -> CSV path
%     .candidate_csv_path       - shared durable homolog CSV in kcatData
%     .candidate_run_csv_path   - run-local shared homolog CSV snapshot
%     .n_candidates_per_enzyme  - table(rxn, nOrganisms) [approx]
%     .qc                       - aggregated QC table (nOrganisms, nProteins,
%                                 nSubstrates, nPredictions per rxn)
%     .phase_logs               - cellstr of structured log messages
%     .elapsed_seconds          - duration

    if nargin < 2 || isempty(predictors)
        predictors = {'UniKP', 'CatPred'};
    end
    if nargin < 3 || isempty(options)
        options = struct();
    end
    if nargin < 4 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.')
        end
    end
    if ~iscell(predictors)
        predictors = cellstr(predictors);
    end
    predictors = canonicalPredictors(predictors);

    % Default injection points ------------------------------------------------
    if ~isfield(options, 'CacheDir') || isempty(options.CacheDir)
        options.CacheDir = fullfile(parameters.reconstructionDir, '.cache', 'uniprot');
    end
    if ~isfield(options, 'dockerFn') || isempty(options.dockerFn)
        options.dockerFn = @defaultDockerFn;
    end
    if ~isfield(options, 'Logger') || isempty(options.Logger)
        options.Logger = @defaultLogger;
    end
    if ~isfield(options, 'organismName') || isempty(options.organismName)
        if isfield(ecModel, 'id') && ~isempty(ecModel.id)
            options.organismName = char(ecModel.id);
        else
            options.organismName = 'ecModel';
        end
        options.organismName = regexprep(options.organismName, '[^\w\-]', '_');
    end
    if ~isfield(options, 'writeDlInput') || isempty(options.writeDlInput)
        options.writeDlInput = @defaultWriteDlInput;
    end
    if ~isfield(options, 'buildSeeds') || isempty(options.buildSeeds)
        options.buildSeeds = @buildOkoReactionSeeds;
    end
    if ~isfield(options, 'fetchEnzymes') || isempty(options.fetchEnzymes)
        options.fetchEnzymes = @fetchCrossSpeciesEnzymes;
    end
    if ~isfield(options, 'buildXlInput') || isempty(options.buildXlInput)
        options.buildXlInput = @buildCrossSpeciesDlInput;
    end
    if ~isfield(options, 'aggregate') || isempty(options.aggregate)
        options.aggregate = @aggregateOkoIntervals;
    end
    if ~isfield(options, 'formatCsv') || isempty(options.formatCsv)
        options.formatCsv = @formatOkoPredictionsCsv;
    end
    if ~isfield(options, 'MaxHomologs') || isempty(options.MaxHomologs)
        options.MaxHomologs = 100;
    end
    if ~isfield(options, 'Resume'), options.Resume = true; end
    if ~isfield(options, 'CheckpointEvery') || isempty(options.CheckpointEvery)
        options.CheckpointEvery = 100;
    end

    t0 = tic;
    sourceDir = fullfile(parameters.reconstructionDir, 'kcatData');
    if ~exist(sourceDir, 'dir'), mkdir(sourceDir); end
    if ~isfield(options, 'OutputDir') || isempty(options.OutputDir)
        runsRoot = fullfile(sourceDir, 'oko_plus_runs');
        if ~exist(runsRoot, 'dir'), mkdir(runsRoot); end
        if options.Resume
            runName = sprintf('%s_h%d_%s', options.organismName, ...
                options.MaxHomologs, lower(strjoin(predictors, '_')));
            runName = regexprep(runName, '[^\w\-]', '_');
            options.OutputDir = fullfile(runsRoot, runName);
        else
            options.OutputDir = tempname(runsRoot);
        end
    end
    outputDir = options.OutputDir;
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end

    predictorPaths = struct();
    predictorInputPaths = struct();
    predictorCheckpointInputPaths = struct();
    predictorOutputPaths = struct();
    predictorCandidatePaths = struct();
    predictorCandidateRunPaths = struct();
    qcCombined = table();
    phaseLogs = {};
    failures = {};

    % Prepare requested model inputs first and select one as the canonical
    % source of reaction/protein/substrate seeds. Predictor inputs differ
    % only in downstream formatting (CatPred adds pdbpath), so homolog
    % discovery must not be repeated per predictor.
    availableInputs = struct();
    sharedSeedCsv = '';
    sharedSeedPredictor = '';
    for p = 1:numel(predictors)
        predictor = predictors{p};
        dlInputCsv = fullfile(sourceDir, [predictor '_input.csv']);
        if ~exist(dlInputCsv, 'file')
            try
                options.writeDlInput(ecModel, predictor, parameters);
            catch err
                options.Logger('writeInputFile failed (%s): %s', predictor, err.message);
                failures{end+1} = sprintf('%s/writeInputFile: %s', predictor, err.message); %#ok<AGROW>
                continue;
            end
        end
        if ~exist(dlInputCsv, 'file')
            options.Logger('DL input CSV missing for %s: %s', predictor, dlInputCsv);
            failures{end+1} = sprintf('%s/input missing', predictor); %#ok<AGROW>
            continue;
        end
        availableInputs.(predictor) = dlInputCsv;
        if isempty(sharedSeedCsv)
            sharedSeedCsv = dlInputCsv;
            sharedSeedPredictor = predictor;
        end
    end
    if isempty(sharedSeedCsv)
        error('buildOkoPlusIntervals:NoInputs', '%s', strjoin(failures, ' | '));
    end

    seeds = options.buildSeeds(sharedSeedCsv, sharedSeedPredictor);
    if isempty(seeds)
        error('buildOkoPlusIntervals:NoSeeds', ...
            'No OKO+ seeds extracted from %s.', sharedSeedCsv);
    end
    options.Logger('Shared seeds: %d reactions', height(seeds));

    durableCandidateCsv = fullfile(sourceDir, 'OKOplus_candidates.csv');
    checkpointMat = fullfile(sourceDir, 'OKOplus_candidates_checkpoint.mat');
    fetchOpts = struct('CacheDir', options.CacheDir, ...
                       'MaxHomologs', options.MaxHomologs, ...
                       'CheckpointCsv', durableCandidateCsv, ...
                       'CheckpointMat', checkpointMat, ...
                       'CheckpointEvery', options.CheckpointEvery, ...
                       'Resume', options.Resume);
    candidates = options.fetchEnzymes(seeds, fetchOpts);
    if isempty(candidates)
        error('buildOkoPlusIntervals:NoCandidates', ...
            'No cross-species candidates were found for the shared seed set.');
    end
    options.Logger('Shared candidates: %d total homolog rows', height(candidates));
    candidateName = 'OKOplus_candidates.csv';
    runCandidateCsv = fullfile(outputDir, candidateName);
    writetable(candidates, runCandidateCsv);
    atomicWriteTable(candidates, durableCandidateCsv);
    options.Logger('Shared candidates saved: %s', durableCandidateCsv);
    candidateSignature = candidatesSignature(candidates);
    for p = 1:numel(predictors)
        predictorCandidatePaths.(predictors{p}) = durableCandidateCsv;
        predictorCandidateRunPaths.(predictors{p}) = runCandidateCsv;
    end

    for p = 1:numel(predictors)
        predictor = predictors{p};
        options.Logger('=== predictor %s ===', predictor);
        if ~isfield(availableInputs, predictor), continue; end

        % 4. Build the durable predictor input once. Candidate checkpoints
        % intentionally do not rebuild these potentially very large CSVs.
        % A successful build gets a signature marker, allowing an interrupted
        % run to reuse it when the exact candidate set is unchanged.
        xlInputCsv = fullfile(sourceDir, [predictor '_OKOplus_input.csv']);
        inputMarker = fullfile(sourceDir, ...
            [predictor '_OKOplus_input_complete.mat']);
        if options.Resume && artifactIsReusable( ...
                inputMarker, xlInputCsv, candidateSignature)
            options.Logger('Predictor input checkpoint reused: %s', predictor);
        elseif options.Resume && predictorInputMatches( ...
                xlInputCsv, seeds, candidates)
            saveArtifactMarker(inputMarker, candidateSignature);
            options.Logger('Existing predictor input verified and adopted: %s', predictor);
        else
            try
                builtInputCsv = options.buildXlInput( ...
                    seeds, candidates, sourceDir, predictor);
                if ~isempty(builtInputCsv)
                    xlInputCsv = builtInputCsv;
                end
                if exist(xlInputCsv, 'file')
                    saveArtifactMarker(inputMarker, candidateSignature);
                end
            catch err
                options.Logger('buildCrossSpeciesDlInput failed (%s): %s', ...
                    predictor, err.message);
                failures{end+1} = sprintf('%s/input: %s', ...
                    predictor, err.message); %#ok<AGROW>
                continue;
            end
        end
        if isempty(xlInputCsv) || ~exist(xlInputCsv, 'file')
            options.Logger('buildCrossSpeciesDlInput produced no CSV for %s', predictor);
            continue;
        end
        predictorCheckpointInputPaths.(predictor) = xlInputCsv;

        % Preserve a human-readable run-local copy as well as the fixed name
        % consumed by each predictor container.
        runOkoInput = fullfile(outputDir, [predictor '_OKOplus_input.csv']);
        stagedInput = fullfile(outputDir, [predictor '_input.csv']);
        stagedMarker = fullfile(outputDir, ...
            [predictor '_input_staged_complete.mat']);
        stagedReusable = options.Resume && ...
            artifactIsReusable(stagedMarker, stagedInput, candidateSignature) && ...
            exist(runOkoInput, 'file');
        if stagedReusable
            options.Logger('Run-local predictor input reused: %s', predictor);
        else
            copyfile(xlInputCsv, runOkoInput);
            % Predictor containers intentionally use fixed filenames. Stage
            % the input under that name without touching reconstruction data.
            copyfile(xlInputCsv, stagedInput);
            saveArtifactMarker(stagedMarker, candidateSignature);
        end
        predictorInputPaths.(predictor) = stagedInput;

        % 5. Docker prediction (or mock). A signature marker allows a rerun
        % to skip a predictor that already completed for this exact shared
        % candidate set.
        predictedCsv = fullfile(outputDir, [predictor '.csv']);
        predictionMarker = fullfile(outputDir, [predictor '_prediction_complete.mat']);
        if options.Resume && artifactIsReusable( ...
                predictionMarker, predictedCsv, candidateSignature)
            options.Logger('Prediction checkpoint reused: %s', predictor);
        else
            try
                options.dockerFn(predictor, outputDir);
                saveArtifactMarker(predictionMarker, candidateSignature);
            catch err
                errMsg = err.message;
                if contains(errMsg, 'docker', 'IgnoreCase', true)
                    options.Logger('Docker missing for %s: %s', predictor, errMsg);
                    failures{end+1} = sprintf('%s/docker: %s', predictor, errMsg); %#ok<AGROW>
                    continue;
                end
                options.Logger('ExecutePrediction failed (%s): %s', predictor, errMsg);
                failures{end+1} = sprintf('%s/predict: %s', predictor, errMsg); %#ok<AGROW>
                continue;
            end
        end

        % 6. Restore original enzymes and aggregate each enzyme-reaction pair
        if ~exist(predictedCsv, 'file')
            options.Logger('Predictor output CSV missing for %s: %s', ...
                predictor, predictedCsv);
            continue;
        end
        predictorOutputPaths.(predictor) = predictedCsv;
        try
            if isequal(options.aggregate, @aggregateOkoIntervals)
                [aggregated, qc] = options.aggregate( ...
                    predictedCsv, seeds, predictor, candidates);
            else
                % Preserve the existing three-argument injection contract.
                [aggregated, qc] = options.aggregate(predictedCsv, seeds, predictor);
            end
        catch err
            options.Logger('aggregateOkoIntervals failed (%s): %s', ...
                predictor, err.message);
            failures{end+1} = sprintf('%s/aggregate: %s', predictor, err.message); %#ok<AGROW>
            continue;
        end

        % 7. Format 6-column CSV
        try
            csvName = sprintf('%s_%s_kcat_preds.csv', ...
                options.organismName, lower(predictor));
            csvPath = fullfile(outputDir, csvName);
            options.formatCsv(aggregated, csvPath);
            predictorPaths.(predictor) = csvPath;
        catch err
            options.Logger('formatOkoPredictionsCsv failed (%s): %s', ...
                predictor, err.message);
            failures{end+1} = sprintf('%s/format: %s', predictor, err.message); %#ok<AGROW>
        end

        % Aggregate QC across predictors
        if ~isempty(qc)
            qcCombined = [qcCombined; qc]; %#ok<AGROW>
        end
    end

    if isempty(fieldnames(predictorPaths))
        if isempty(failures), failures = {'No predictor produced an interval CSV.'}; end
        error('buildOkoPlusIntervals:NoArtifacts', '%s', strjoin(failures, ' | '));
    end

    % nCandidates summary (rxn, nOrganisms)
    nCandidatesTable = qcCombined;
    if ~isempty(nCandidatesTable) && any(strcmp(nCandidatesTable.Properties.VariableNames, 'nOrganisms'))
        keep = {'rxn'};
        if ismember('uniprot', nCandidatesTable.Properties.VariableNames)
            keep{end+1} = 'uniprot'; %#ok<AGROW>
        end
        keep{end+1} = 'nOrganisms';
        nCandidatesTable = nCandidatesTable(:, keep);
        nCandidatesTable.Properties.VariableNames{'nOrganisms'} = 'nHomologs';
    end

    result = struct( ...
        'predictor_csv_paths', predictorPaths, ...
        'predictor_input_paths', predictorInputPaths, ...
        'predictor_checkpoint_input_paths', predictorCheckpointInputPaths, ...
        'predictor_output_paths', predictorOutputPaths, ...
        'predictor_candidate_paths', predictorCandidatePaths, ...
        'predictor_candidate_run_paths', predictorCandidateRunPaths, ...
        'candidate_csv_path', durableCandidateCsv, ...
        'candidate_checkpoint_path', checkpointMat, ...
        'candidate_run_csv_path', runCandidateCsv, ...
        'n_candidates_per_enzyme', nCandidatesTable, ...
        'qc', qcCombined, ...
        'phase_logs', {phaseLogs}, ...
        'run_dir', outputDir, ...
        'elapsed_seconds', toc(t0));
end

function signature = candidatesSignature(candidates)
    parts = string(candidates.ReactionName) + "|" + ...
            string(candidates.OriginalProteinID) + "|" + ...
            string(candidates.accession);
    payload = unicode2native(char(strjoin(parts, newline)), 'UTF-8');
    md = java.security.MessageDigest.getInstance('SHA-256');
    hash = md.digest(uint8(payload));
    signature = lower(reshape(dec2hex(hash), 1, []));
end

function reusable = artifactIsReusable(markerPath, outputPath, signature)
    reusable = false;
    if ~exist(markerPath, 'file') || ~exist(outputPath, 'file'), return; end
    try
        state = load(markerPath, 'candidateSignature');
        reusable = isfield(state, 'candidateSignature') && ...
            strcmp(state.candidateSignature, signature);
    catch
        reusable = false;
    end
end

function saveArtifactMarker(markerPath, candidateSignature)
    tempPath = [tempname(fileparts(markerPath)) '.mat'];
    cleanup = onCleanup(@()deleteIfPresent(tempPath)); %#ok<NASGU>
    save(tempPath, 'candidateSignature', '-v7');
    [ok, message] = movefile(tempPath, markerPath, 'f');
    if ~ok, error('runOkoPipeline:MarkerSaveFailed', '%s', message); end
end

function matches = predictorInputMatches(csvPath, seeds, candidates)
%PREDICTORINPUTMATCHES Safely adopt large inputs made before markers existed.
% Compare the exact multiplicity of each (reaction, homolog accession) pair;
% substrate expansion determines that multiplicity for every seed.
    matches = false;
    if ~exist(csvPath, 'file'), return; end
    try
        expectedKeys = strings(0,1);
        for i = 1:height(seeds)
            rxn = string(seeds.ReactionName{i});
            orig = string(seeds.OriginalProteinID{i});
            mask = string(candidates.ReactionName) == rxn;
            if ismember('OriginalProteinID', candidates.Properties.VariableNames)
                mask = mask & string(candidates.OriginalProteinID) == orig;
            end
            accessions = string(candidates.accession(mask));
            nSubstrates = numel(seeds.Substrates{i});
            expectedKeys = [expectedKeys; ...
                repelem(rxn + "|" + accessions, nSubstrates)]; %#ok<AGROW>
        end

        importOpts = detectImportOptions(csvPath, 'VariableNamingRule', 'preserve');
        required = {'ReactionName','ProteinID'};
        if ~all(ismember(required, importOpts.VariableNames)), return; end
        importOpts.SelectedVariableNames = required;
        actual = readtable(csvPath, importOpts);
        actualKeys = string(actual.ReactionName) + "|" + string(actual.ProteinID);
        if numel(actualKeys) ~= numel(expectedKeys), return; end
        matches = isequal(sort(actualKeys), sort(expectedKeys));
    catch
        matches = false;
    end
end

function predictors = canonicalPredictors(predictors)
    allowed = {'DLKcat','UniKP','CatPred'};
    canonical = cell(size(predictors));
    for i = 1:numel(predictors)
        hit = find(strcmpi(char(predictors{i}), allowed), 1);
        if isempty(hit)
            error('buildOkoPlusIntervals:InvalidPredictor', ...
                'Predictor must be DLKcat, UniKP or CatPred: %s', ...
                char(predictors{i}));
        end
        canonical{i} = allowed{hit};
    end
    predictors = unique(canonical, 'stable');
end

function atomicWriteTable(T, targetPath)
%ATOMICWRITETABLE Avoid leaving a truncated durable artifact on interruption.
    tempPath = [tempname(fileparts(targetPath)) '.csv'];
    cleanup = onCleanup(@()deleteIfPresent(tempPath)); %#ok<NASGU>
    writetable(T, tempPath);
    [ok, message] = movefile(tempPath, targetPath, 'f');
    if ~ok
        error('runOkoPipeline:CandidateSaveFailed', '%s', message);
    end
end

function deleteIfPresent(path)
    if exist(path, 'file'), delete(path); end
end

% ---------------------------------------------------------------------------
%  Defaults
% ---------------------------------------------------------------------------

function defaultWriteDlInput(ecModel, predictor, params)
%DEFAULTWRITEDLINPUT Generate the standard DL input CSV via writeInputFile.
    writeInputFile(ecModel, predictor, params);
end

function defaultDockerFn(predictor, outputDir)
%DEFAULTDOCKERFN Run the real Docker predictor via ExecutePrediction.
    ExecutePrediction({predictor}, outputDir);
end

function defaultLogger(varargin)
%DEFAULTLOGGER Lightweight fprintf-based logger.
    prefix = '[runOkoPipeline] ';
    if numel(varargin) == 1
        fprintf([prefix '%s\n'], varargin{1});
    else
        fmt = ['%s' varargin{1} '\n'];
        args = [{prefix} varargin(2:end)];
        fprintf(fmt, args{:});
    end
end
