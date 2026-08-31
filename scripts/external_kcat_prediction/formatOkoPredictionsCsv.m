function csvPath = formatOkoPredictionsCsv(aggregatedTable, outputPath)
%FORMATOKOPREDICTIONSCSV Write OKO-compatible CSV (6 columns, legacy format).
%
%   The output strictly matches the legacy OKO+ solver's expected schema
%   from OKO-main/Ecoli/ecoli_kcat_preds.csv:
%       rxn, uniprot, min, max, mean, mode
%
%   IMPORTANT: mean and mode are compatibility metadata only and do NOT
%   participate in OKO+ MILP bounds. The solver reads only min/max.
%
%   csvPath = formatOkoPredictionsCsv(aggregatedTable, outputPath)
%
%   See: docs/superpowers/specs/2026-08-09-oko-plus-homolog-interval-pipeline-design.md
%        section 4.4

    requiredCols = {'rxn', 'uniprot', 'min', 'max', 'mean', 'mode'};
    if ~all(ismember(requiredCols, aggregatedTable.Properties.VariableNames))
        error('formatOkoPredictionsCsv:missingColumns', ...
              'aggregatedTable must contain columns: %s', strjoin(requiredCols, ', '));
    end

    out = aggregatedTable(:, requiredCols);

    % Ensure outputPath has .csv extension
    [outDir, baseName, ext] = fileparts(outputPath);
    if ~strcmp(ext, '.csv')
        outputPath = fullfile(outDir, [baseName '.csv']);
    end
    if ~isempty(outDir) && ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    writetable(out, outputPath);
    csvPath = outputPath;
end