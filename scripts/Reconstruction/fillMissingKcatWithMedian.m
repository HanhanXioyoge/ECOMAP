function [model, T] = fillMissingKcatWithMedian(model, maxRows, fillMedian)
    % fillMissingKcatWithMedianECOMAP
    % List reactions in model.enzymeConstraints with kcat == 0 and show,
    % for each reaction, the enzyme indices involved plus the paired
    % enzyme IDs and gene names at those indices. Additionally, fill
    % missing kcat values with the median kcat if specified.
    %
    % Output:
    %   model        - Updated ECOMAP model with filled kcat values
    %   T            - Table with columns:
    %      rxnIdxEC    - index into model.enzymeConstraints.rxns
    %      rxnID_EC    - reaction ID in enzymeConstraints
    %      rxnName     - model.rxnNames (if resolvable)
    %      nEnz        - number of enzymes linked to this reaction
    %      pairs       - cell array of strings "enzymes{j} <-> genes{j} (idx j)"
    %
    % Inputs:
    %   model        - ECOMAP model with enzyme constraints (contains enzymeConstraints structure)
    %   maxRows      - Maximum number of rows to display (optional)
    %   fillMedian   - Boolean flag indicating whether to fill missing kcat with median value

    if nargin < 3, fillMedian = false; end  % Default: do not fill with median if not specified
    if nargin < 2 || isempty(maxRows), maxRows = Inf; end

    ec = model.enzymeConstraints;
    nEC = numel(ec.rxns);

    % Ensure kcat vector exists
    if ~isfield(ec, 'kcat') || numel(ec.kcat) ~= nEC
        ec.kcat = zeros(nEC, 1);
    end

    % Basic field checks
    req = {'rxns', 'kcat', 'enzymes', 'genes', 'rxnEnzMat'};
    for k = 1:numel(req)
        if ~isfield(ec, req{k})
            error('model.enzymeConstraints.%s is missing.', req{k});
        end
    end

    % Find reactions with kcat == 0
    missIdx = find(ec.kcat == 0);
    if isempty(missIdx)
        T = table();
        fprintf('[fillMissingKcatWithMedianECOMAP] All kcat values are filled (no zeros found).\n');
        return;
    end

    % Fill missing kcat values with median if requested
    if fillMedian
        medianKcat = median(ec.kcat(ec.kcat > 0));  % Calculate median of non-zero kcat values
        if isempty(medianKcat)
            warning('No non-zero kcat values available to calculate median.');
        else
            ec.kcat(missIdx) = medianKcat;
            model.enzymeConstraints.kcat = ec.kcat;
            model.enzymeConstraints.source(missIdx) = {'median'};  % Mark source as "median fill"
            fprintf('[fillMissingKcatWithMedianECOMAP] Missing kcat values filled with median: %f\n', medianKcat);
        end
    end

    % Map EC rxns back to model.rxnNames (best effort)
    rxnNames = repmat({''}, numel(missIdx), 1);
    try
        ecIDs = ec.rxns(missIdx);
        switch getfield(ec, 'ecModeltype') %#ok<GFLD>
            case 'basic'
                baseIDs = extractAfter(ecIDs, 4); % strip '001_' style prefixes
            otherwise
                baseIDs = ecIDs;
        end
        [ok, idxModel] = ismember(baseIDs, model.rxns);
        rxnNames(ok) = model.rxnNames(idxModel(ok));
    catch
        % silently ignore if model.rxnNames is not available
    end

    % Build pairs "enzyme <-> gene (idx)"
    pairsCol = cell(numel(missIdx), 1);
    nEnzCol  = zeros(numel(missIdx), 1);

    for i = 1:numel(missIdx)
        r = missIdx(i);
        enzMask = ec.rxnEnzMat(r, :) ~= 0;
        enzIdx  = find(enzMask);
        nEnzCol(i) = numel(enzIdx);

        if isempty(enzIdx)
            pairsCol{i} = {'<no enzymes linked>'};
        else
            % Ensure index-aligned access
            enzIDs  = ec.enzymes(enzIdx);
            geneIDs = ec.genes(enzIdx);
            line    = arrayfun(@(j) sprintf('%s <-> %s (idx %d)', ...
                                enzIDs{j}, geneIDs{j}, enzIdx(j)), ...
                                (1:numel(enzIdx))', 'UniformOutput', false);
            pairsCol{i} = line;
        end
    end

    % Assemble table and optionally truncate
    T = table(missIdx, ec.rxns(missIdx), rxnNames, nEnzCol, pairsCol, ...
        'VariableNames', {'rxnIdxEC','rxnID_EC','rxnName','nEnz','pairs'});

    if isfinite(maxRows) && height(T) > maxRows
        T = T(1:maxRows, :);
    end

    % Pretty console preview (first few rows)
    previewN = min(5, height(T));
    fprintf('[fillMissingKcatWithMedianECOMAP] %d reactions with kcat==0. Showing first %d:\n', ...
            numel(missIdx), previewN);
    for i = 1:previewN
        fprintf('  #%d  %s  | %s  | nEnz=%d\n', ...
            T.rxnIdxEC(i), T.rxnID_EC{i}, T.rxnName{i}, T.nEnz(i));
        p = T.pairs{i};
        for j = 1:min(numel(p), 3)
            fprintf('      - %s\n', p{j});
        end
        if numel(p) > 3, fprintf('      ... (%d more)\n', numel(p)-3); end
    end
end
