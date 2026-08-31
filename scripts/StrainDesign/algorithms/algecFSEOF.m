function fseof = algecFSEOF(model, prodTargetRxn, csRxn, nSteps, outputFile, filePath, parameters)
%ALGECFSEOF Run ecFSEOF for an ECOMAP enzyme-constrained model.
%
%   fseof = algecFSEOF(model, prodTargetRxn, csRxn, nSteps, ...
%                      outputFile, filePath, parameters)
%
% Inputs
%   model          ECOMAP ecGEM. model.type must identify ECOMAP and
%                  model.enzymeConstraints.ecModeltype must be basic,
%                  isozyme, or integrated.
%   prodTargetRxn  Product reaction ID; a secretion-positive exchange
%                  reaction is recommended.
%   csRxn          Main carbon-source uptake reaction ID. Its bound is
%                  checked but is not changed by this function.
%   nSteps         Number of enforced-production steps (default 16; >= 2).
%   outputFile     Logical flag controlling file output (default false).
%   filePath       Base output directory (default parameters.designDir).
%   parameters     Project parameter struct containing parameters.bioRxn.
%
% Files are written directly into filePath. Their names identify the
% algorithm, ecModel structure, product, and result type, for example:
%   ecFSEOF_integrated_EX_ala__L_e_genes.tsv
%   ecFSEOF_integrated_EX_ala__L_e_rxns.tsv
%   ecFSEOF_integrated_EX_ala__L_e_transporter.tsv
%   ecFSEOF_integrated_EX_ala__L_e_result.mat
%
% Output fields
%   modelID, modelType, prodTargetRxn, carbonSourceRxn, outputDir,
%   alpha, v_matrix, rxnTargets, transportTargets, and geneTargets.

    if nargin < 7 || isempty(parameters)
        parameters = ParameterManager.getParams();
    end
    if ~isstruct(parameters) || ~isfield(parameters, 'bioRxn') || ...
            isempty(parameters.bioRxn)
        error('algecFSEOF:MissingBiomassReaction', ...
            'parameters.bioRxn is required.');
    end
    if nargin < 4 || isempty(nSteps)
        nSteps = 16;
    end
    if nargin < 5 || isempty(outputFile)
        outputFile = false;
    end
    if nargin < 6 || isempty(filePath)
        if ~isfield(parameters, 'designDir') || isempty(parameters.designDir)
            error('algecFSEOF:MissingDesignDirectory', ...
                'parameters.designDir is required when filePath is omitted.');
        end
        filePath = parameters.designDir;
    end

    prodTargetRxn = char(string(prodTargetRxn));
    csRxn = char(string(csRxn));
    bioRxn = char(string(parameters.bioRxn));
    validateattributes(nSteps, {'numeric'}, ...
        {'scalar', 'integer', '>=', 2, 'finite'}, mfilename, 'nSteps');
    validateattributes(outputFile, {'logical', 'numeric'}, ...
        {'scalar'}, mfilename, 'outputFile');
    outputFile = logical(outputFile);

    [modelID, ecModeltype] = validateEcomapEcModel(model);
    prodTargetRxnIdx = resolveReaction(model, prodTargetRxn, 'product target');
    csRxnIdx = resolveReaction(model, csRxn, 'carbon source');
    bioRxnIdx = resolveReaction(model, bioRxn, 'biomass');

    outputDir = '';
    outputPrefix = '';
    if outputFile
        if ~(ischar(filePath) || (isstring(filePath) && isscalar(filePath)))
            error('algecFSEOF:BadOutputDirectory', ...
                'filePath must be a character vector or scalar string directory.');
        end
        outputDir = char(filePath);
        outputPrefix = resultFilePrefix(ecModeltype, prodTargetRxn);
        if ~exist(outputDir, 'dir')
            [ok, message] = mkdir(outputDir);
            if ~ok
                error('algecFSEOF:CreateOutputDirectory', ...
                    'Could not create output directory "%s": %s', ...
                    outputDir, message);
            end
        end
    end

    % Standardize gene rules before reaction/gene target mapping.
    [grRules, rxnGeneMat] = standardizeGrRules(model, true);
    model.grRules = grRules;
    model.rxnGeneMat = rxnGeneMat;

    % Baseline: maximize growth. csRxn is inspected only; its bounds remain
    % exactly as supplied by the caller.
    growthModel = setParam(model, 'obj', bioRxn, 1);
    growthSol = solveLP(growthModel);
    assertSolution(growthSol, numel(model.rxns), 'baseline growth maximization');
    if model.lb(csRxnIdx) < growthSol.x(csRxnIdx)
        printOrange(['WARNING: Carbon source lower bound was set to ' ...
            num2str(model.lb(csRxnIdx)) ', but the uptake rate after ' ...
            'model optimization is ' num2str(growthSol.x(csRxnIdx)) '.\n']);
    end
    iniTarget = growthSol.x(prodTargetRxnIdx);

    % Product maximum and enforced-production scan.
    productModel = setParam(model, 'obj', prodTargetRxn, 1);
    productSol = solveLP(productModel, 1);
    assertSolution(productSol, numel(model.rxns), 'product maximization');
    theoreticalTargetMax = productSol.x(prodTargetRxnIdx);
    if ~isfinite(theoreticalTargetMax) || theoreticalTargetMax <= 0
        error('algecFSEOF:NonPositiveTargetMaximum', ...
            ['The maximum flux of product reaction "%s" is %g. ', ...
             'Check the reaction direction and model constraints.'], ...
            prodTargetRxn, theoreticalTargetMax);
    end
    maxTarget = theoreticalTargetMax * 0.9;
    fluxTolerance = 1e-9;
    if maxTarget <= iniTarget + fluxTolerance
        error('algecFSEOF:InvalidScanRange', ...
            ['The 90%% product maximum (%g) must exceed the product flux ', ...
             'at maximum growth (%g).'], maxTarget, iniTarget);
    end

    alpha = linspace(iniTarget, maxTarget, nSteps);
    vMatrixAll = zeros(numel(model.rxns), nSteps);
    scanModel = model;
    progressbar('Flux Scanning with Enforced Objective Function')
    for i = 1:nSteps
        scanModel = setParam(scanModel, 'eq', prodTargetRxn, alpha(i));
        scanModel.lb(bioRxnIdx) = 0;
        scanModel = setParam(scanModel, 'obj', bioRxn, 1);
        sol = solveLP(scanModel, 1);
        assertSolution(sol, numel(model.rxns), ...
            sprintf('ecFSEOF step %d/%d (product flux %g)', ...
                    i, nSteps, alpha(i)));
        vMatrixAll(:, i) = sol.x(:);
        progressbar(i / nSteps)
    end
    progressbar(1)

    % Keep gene-associated reactions, excluding synthetic standard genes.
    withGR = ~cellfun(@isempty, model.grRules);
    withGR(contains(model.grRules, 'standard')) = false;
    targetRxns = model.rxns(withGR);
    vMatrix = vMatrixAll(withGR, :);
    rxnGeneM = model.rxnGeneMat(withGR, :);

    nonzeroFlux = any(abs(vMatrix) > fluxTolerance, 2);
    targetRxns = targetRxns(nonzeroFlux);
    vMatrix = vMatrix(nonzeroFlux, :);
    rxnGeneM = rxnGeneM(nonzeroFlux, :);

    % Identify reactions whose absolute flux changes monotonically.
    nCandidates = numel(targetRxns);
    slopeRxns = zeros(nCandidates, 1);
    targets = false(nCandidates, 1);
    targetType = cell(nCandidates, 1);
    scanRange = alpha(end) - alpha(1);
    for i = 1:nCandidates
        fluxDelta = diff(abs(vMatrix(i, :)));
        if all(fluxDelta > fluxTolerance)
            targets(i) = true;
            slopeRxns(i) = abs(vMatrix(i, end) - vMatrix(i, 1)) / scanRange;
            targetType{i} = 'OE';
        elseif all(fluxDelta < -fluxTolerance)
            targets(i) = true;
            slopeRxns(i) = abs(vMatrix(i, end) - vMatrix(i, 1)) / scanRange;
            if abs(vMatrix(i, end)) <= fluxTolerance
                targetType{i} = 'KO';
            else
                targetType{i} = 'KD';
            end
        end
    end

    targetRxns = targetRxns(targets);
    vMatrix = vMatrix(targets, :);
    rxnGeneM = rxnGeneM(targets, :);
    slopeRxns = slopeRxns(targets);
    targetType = targetType(targets);

    [slopeRxns, order] = sort(slopeRxns, 'descend');
    targetRxns = targetRxns(order);
    vMatrix = vMatrix(order, :);
    rxnGeneM = rxnGeneM(order, :);
    targetType = targetType(order);

    % Preserve the original ecFSEOF policy: retain the upper slope quartile.
    if ~isempty(slopeRxns)
        keep = slopeRxns > quantile(slopeRxns, 0.75);
        targetRxns = targetRxns(keep);
        vMatrix = vMatrix(keep, :);
        rxnGeneM = rxnGeneM(keep, :);
        slopeRxns = slopeRxns(keep);
        targetType = targetType(keep);
    end

    % Build gene-level targets.
    geneMask = sum(rxnGeneM, 1) > 0;
    genes = model.genes(geneMask);
    rxnGeneM = rxnGeneM(:, geneMask);
    slopeGenes = zeros(numel(genes), 1);
    targetTypeGenes = cell(numel(genes), 1);
    essentiality = repmat({'nonessential'}, numel(genes), 1);
    essentialityModel = setParam(model, 'obj', bioRxn, 1);

    progressbar('Checking for gene essentiality')
    for i = 1:numel(genes)
        tempModel = geneKnockoutModel(essentialityModel, genes{i}, ecModeltype);
        solKO = solveLP(tempModel);
        isEssential = ~hasValidSolution(solKO, numel(model.rxns)) || ...
            solKO.x(bioRxnIdx) < 1e-8;
        if isEssential
            essentiality{i} = 'essential';
        end

        rxnsForGene = find(rxnGeneM(:, i) > 0);
        actions = unique(targetType(rxnsForGene));
        targetTypeGenes{i} = reconcileGeneAction(actions, isEssential);
        slopeGenes(i) = mean(slopeRxns(rxnsForGene));
        progressbar(i / max(1, numel(genes)))
    end
    progressbar(1)

    [slopeGenes, order] = sort(slopeGenes, 'descend');
    genes = genes(order);
    targetTypeGenes = targetTypeGenes(order);
    essentiality = essentiality(order);

    % Metadata makes in-memory results self-identifying as well.
    fseof.modelID = modelID;
    fseof.modelType = ecModeltype;
    fseof.prodTargetRxn = prodTargetRxn;
    fseof.carbonSourceRxn = csRxn;
    fseof.biomassRxn = bioRxn;
    fseof.outputDir = outputDir;
    fseof.outputPrefix = outputPrefix;
    fseof.alpha = alpha;

    % Exclude integrated protein-usage reactions from reaction reports.
    toKeep = ~startsWith(targetRxns, 'usage_prot_') & ...
        ~startsWith(targetRxns, 'draw_prot_');
    reportedRxns = targetRxns(toKeep);
    reportedFlux = vMatrix(toKeep, :);
    reportedSlopes = slopeRxns(toKeep);
    if isempty(reportedRxns)
        rxnIdx = zeros(0, 1);
        fseof.v_matrix = array2table(zeros(0, nSteps), ...
            'VariableNames', matlab.lang.makeUniqueStrings( ...
                matlab.lang.makeValidName(cellstr(string(alpha)))));
        fseof.rxnTargets = cell(0, 5);
    else
        rxnIdx = getIndexes(model, reportedRxns, 'rxns');
        variableNames = matlab.lang.makeUniqueStrings( ...
            matlab.lang.makeValidName(cellstr(string(alpha))));
        fseof.v_matrix = array2table(reportedFlux, ...
            'VariableNames', variableNames, 'RowNames', model.rxns(rxnIdx));
        fseof.rxnTargets = [model.rxns(rxnIdx), model.rxnNames(rxnIdx), ...
            num2cell(reportedSlopes), model.grRules(rxnIdx), ...
            constructEquations(model, rxnIdx)];
    end

    transportRxns = false(numel(rxnIdx), 1);
    for i = 1:numel(rxnIdx)
        mets = model.metNames(model.S(:, rxnIdx(i)) ~= 0);
        mets(startsWith(mets, 'prot_')) = [];
        transportRxns(i) = numel(mets) ~= numel(unique(mets));
    end
    fseof.transportTargets = fseof.rxnTargets(transportRxns, :);
    fseof.rxnTargets(transportRxns, :) = [];

    if isempty(genes)
        fseof.geneTargets = cell(0, 5);
    else
        [foundGenes, geneIdx] = ismember(lower(string(genes)), ...
            lower(string(model.genes)));
        if ~all(foundGenes)
            error('algecFSEOF:GeneMappingFailed', ...
                'Could not map all target genes back to model.genes.');
        end
        shortNames = model.genes(geneIdx);
        if isfield(model, 'geneShortNames') && ...
                numel(model.geneShortNames) == numel(model.genes)
            shortNames = model.geneShortNames(geneIdx);
        end
        fseof.geneTargets = [genes(:), shortNames(:), ...
            num2cell(slopeGenes), targetTypeGenes(:), essentiality(:)];
    end

    if outputFile
        writeOutputs(fseof, outputDir, outputPrefix);
        fprintf('[algecFSEOF] model=%s, type=%s, target=%s\n', ...
            modelID, ecModeltype, prodTargetRxn);
        fprintf('[algecFSEOF] Results saved to %s\n', outputDir);
    end
end

function [modelID, ecModeltype] = validateEcomapEcModel(model)
    if ~isstruct(model) || ~isscalar(model)
        error('algecFSEOF:InvalidModel', 'model must be a scalar structure.');
    end
    required = {'S', 'rxns', 'lb', 'ub', 'type', 'enzymeConstraints'};
    for i = 1:numel(required)
        if ~isfield(model, required{i}) || isempty(model.(required{i}))
            error('algecFSEOF:NotEcomapEcModel', ...
                'The input is not an ECOMAP ecGEM: missing model.%s.', required{i});
        end
    end
    if size(model.S, 2) ~= numel(model.rxns) || ...
            numel(model.lb) ~= numel(model.rxns) || ...
            numel(model.ub) ~= numel(model.rxns)
        error('algecFSEOF:InvalidModelDimensions', ...
            'model.S, model.rxns, model.lb, and model.ub dimensions are inconsistent.');
    end

    declaredType = lower(strtrim(char(string(model.type))));
    if ~startsWith(declaredType, 'ecomap')
        error('algecFSEOF:NotEcomapEcModel', ...
            'algecFSEOF requires an ECOMAP ecGEM; model.type is "%s".', ...
            char(string(model.type)));
    end
    ec = model.enzymeConstraints;
    if ~isstruct(ec) || ~isfield(ec, 'ecModeltype') || isempty(ec.ecModeltype)
        error('algecFSEOF:MissingEcModelType', ...
            'model.enzymeConstraints.ecModeltype is required.');
    end
    ecModeltype = lower(strtrim(char(string(ec.ecModeltype))));
    supported = {'basic', 'isozyme', 'integrated'};
    if ~ismember(ecModeltype, supported)
        error('algecFSEOF:UnsupportedEcModelType', ...
            'Unsupported ecModeltype "%s". Supported types: %s.', ...
            ecModeltype, strjoin(supported, ', '));
    end
    typeSuffix = regexprep(declaredType, '^ecomap[-_ ]*', '');
    if ~isempty(typeSuffix) && ~strcmp(typeSuffix, ecModeltype)
        error('algecFSEOF:EcModelTypeMismatch', ...
            ['model.type "%s" conflicts with ', ...
             'model.enzymeConstraints.ecModeltype "%s".'], ...
            char(string(model.type)), ecModeltype);
    end

    modelID = '';
    for field = {'id', 'name', 'modelID', 'modelName'}
        if isfield(model, field{1}) && ~isempty(model.(field{1}))
            value = model.(field{1});
            if ischar(value) || (isstring(value) && isscalar(value))
                modelID = char(value);
                break;
            end
        end
    end
    if isempty(strtrim(modelID))
        error('algecFSEOF:MissingModelIdentity', ...
            ['model.id, model.name, model.modelID, or model.modelName is ', ...
             'required to keep output from different models separate.']);
    end
end

function idx = resolveReaction(model, rxnID, role)
    idx = find(strcmp(model.rxns, rxnID));
    if isempty(idx)
        error('algecFSEOF:ReactionNotFound', ...
            'The %s reaction "%s" is not present in model.rxns.', role, rxnID);
    end
    if numel(idx) > 1
        error('algecFSEOF:ReactionNotUnique', ...
            'The %s reaction "%s" occurs %d times.', role, rxnID, numel(idx));
    end
end

function assertSolution(sol, nRxns, context)
    if ~hasValidSolution(sol, nRxns)
        error('algecFSEOF:InvalidSolution', ...
            'No valid LP solution was returned for %s.', context);
    end
end

function valid = hasValidSolution(sol, nRxns)
    valid = isstruct(sol) && isfield(sol, 'x') && ...
        numel(sol.x) == nRxns && all(isfinite(sol.x));
    if valid && isfield(sol, 'stat')
        valid = ~isequal(sol.stat, -1);
    end
end

function tempModel = geneKnockoutModel(model, gene, ecModeltype)
    if strcmp(ecModeltype, 'integrated')
        ec = model.enzymeConstraints;
        required = {'genes', 'enzymes'};
        for i = 1:numel(required)
            if ~isfield(ec, required{i})
                error('algecFSEOF:MissingEnzymeMapping', ...
                    'Integrated model is missing enzymeConstraints.%s.', ...
                    required{i});
            end
        end
        enzymeMask = strcmpi(ec.genes, gene);
        usageRxns = strcat('usage_prot_', ec.enzymes(enzymeMask));
        usageRxns = usageRxns(ismember(usageRxns, model.rxns));
        if ~isempty(usageRxns)
            tempModel = setParam(model, 'eq', usageRxns, 0);
            return;
        end
    end

    % basic/isozyme models have no usage_prot_* reactions. Integrated
    % models can also contain GPR genes that were not retained in the
    % enzymeConstraints mapping, so use the same GPR-aware fallback.
    if exist('removeGenes', 'file') ~= 2
        error('algecFSEOF:MissingGeneDeletionFunction', ...
            ['Gene "%s" has no usable protein-usage reaction. ', ...
             'removeGenes from RAVEN is required.'], gene);
    end
    tempModel = removeGenes(model, {gene}, false, false, false);
end

function action = reconcileGeneAction(actions, isEssential)
    if any(strcmp(actions, 'OE'))
        action = 'OE';
    elseif isEssential
        action = 'KD';
    elseif any(strcmp(actions, 'KO'))
        action = 'KO';
    else
        action = 'KD';
    end
end

function prefix = resultFilePrefix(ecModeltype, targetRxn)
    prefix = ['ecFSEOF_' safePathComponent(ecModeltype) '_' ...
              safePathComponent(targetRxn)];
end

function safe = safePathComponent(value)
    value = char(string(value));
    chunks = cell(1, numel(value));
    for i = 1:numel(value)
        ch = value(i);
        if ~isempty(regexp(ch, '[A-Za-z0-9_-]', 'once'))
            chunks{i} = ch;
        else
            chunks{i} = sprintf('~%04X~', double(ch));
        end
    end
    safe = [chunks{:}];
    if isempty(safe)
        safe = 'unnamed';
    end
end

function writeOutputs(fseof, outputDir, outputPrefix)
    genesTable = cell2table(fseof.geneTargets, ...
        'VariableNames', {'gene_IDs', 'gene_names', 'slope', ...
                          'action', 'essentiality'});
    rxnsTable = cell2table(fseof.rxnTargets, ...
        'VariableNames', {'rxn_IDs', 'rxnNames', 'slope', ...
                          'grRules', 'rxn_formula'});
    transportTable = cell2table(fseof.transportTargets, ...
        'VariableNames', {'rxn_IDs', 'rxnNames', 'slope', ...
                          'grRules', 'rxn_formula'});

    writetable(genesTable, fullfile(outputDir, ...
        [outputPrefix '_genes.tsv']), ...
        'FileType', 'text', 'Delimiter', '\t', 'QuoteStrings', false);
    writetable(rxnsTable, fullfile(outputDir, ...
        [outputPrefix '_rxns.tsv']), ...
        'FileType', 'text', 'Delimiter', '\t', 'QuoteStrings', false);
    writetable(transportTable, ...
        fullfile(outputDir, [outputPrefix '_transporter.tsv']), ...
        'FileType', 'text', 'Delimiter', '\t', 'QuoteStrings', false);
    save(fullfile(outputDir, [outputPrefix '_result.mat']), 'fseof');
end
