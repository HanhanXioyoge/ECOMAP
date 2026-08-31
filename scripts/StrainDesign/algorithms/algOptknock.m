function result = algOptknock(model, biomassRxn, targetRxn, selectedRxnList, numDel, maxSolutions, values, sense, verbFlag, outputFile,parameters)
% algOptknock  Find gene knockout sets using COBRA OptKnock.
%   result = algOptknock(model, biomassRxn, targetRxn, selectedRxnList, ...
%                        numDel, maxSolutions, values, sense, verbFlag, outputFile)
%
%   Wraps COBRA Toolbox OptKnock to find reaction knockout sets that
%   maximize target production while maintaining minimum growth.
%   All outputs are KO targets (knockout targets).
%
%   Re-run overwrites previous optknock.* files in the workspace subfolder.
%   Per-solution knockout reactions are grouped together in the CSV;
%   solutions are visually separated by a blank line.
%
%   Inputs:
%     model            struct, COBRA-compatible model
%     biomassRxn       char, biomass reaction ID
%     targetRxn        char, target reaction ID
%     selectedRxnList  cellstr / numeric, candidates for knockout
%     numDel           int, max number of deletions per solution
%     maxSolutions     int, max number of solutions to enumerate
%     values           double, fraction of WT biomass (default 0.1)
%     sense            char, constraint sense (default 'GE')
%     verbFlag         logical, verbose Command Window logging (default false)
%     outputFile       char, output CSV path (empty -> stdout;
%                                       default -> workspace CS.V path)
%
%   Output:
%     result - struct (algFseof-style 7 fields):
%         config          - actual params used
%         biomassRxn      - biomass reaction ID
%         targetRxn       - target reaction ID
%         outputFile      - CSV output path (or '' if stdout)
%         matFile         - path to optknock_result.mat
%         rows            - struct array, one row per KO hit
%                           (grouped by solutionID, blank-line separator in CSV)
%         targets         - flat unique KO list + appearanceCount + diagnostics
    % ------------------- Parameters & inputs -------------------
    if nargin < 11 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end
    % --- Type squash ---
    biomassRxn = char(biomassRxn);
    targetRxn  = char(targetRxn);

    % --- Defaults (nargin chain, mirroring algFseof) ---
    if nargin < 5
        numDel       = 5;
        maxSolutions = 200;
        values       = 0.1;
        sense        = 'G';
        verbFlag     = false;
    end
    if nargin < 6
        maxSolutions = 200;
    end
    if nargin < 7
        values = 0.1;
    end
    if nargin < 8
        sense = 'G';
    end
    if nargin < 9
        verbFlag = false;
    end
    if nargin < 10
        outputFile = parameters.designDir;
    end
    % --- outputFile default (mirror algFseof) ---
    safeTarget = safeTargetName(targetRxn);
    csvName    = ['optknock_' safeTarget '.csv'];
    matName    = ['optknock_' safeTarget '_result.mat'];
    outputFile = fullfile(outputFile, csvName);
    algoDir    = fileparts(outputFile);
    matFile    = fullfile(algoDir, matName);

    output = ~isempty(outputFile);

    % --- Re-run cleanup (mirror algFseof) ---
    for oldName = {csvName, matName}
        p = fullfile(algoDir, oldName{1});
        if exist(p, 'file'), delete(p); end
    end

    % --- Solver setup (preserved from original) ---
    % changeCobraSolver('gurobi','all');
    model=setParam(model,'obj',biomassRxn,1);
    sol = solveLP(model);
    wtBiomass = sol.f;

    % --- Target max flux (for diagnostics only) ---
    targetMax = NaN;
    try
        m2 = setParam(model, 'obj', targetRxn, 1);
        s2 = solveLP(m2, 1);
        targetMax = s2.f;
    catch
    end

    % --- COBRA's internal solution file path (independent of our .mat) ---
    solutionFileNameTmp = fullfile(algoDir, ...
        ['optknock_' safeTarget '_cobra.sol']);

    % --- Build result struct skeleton (algFseof 7-field shape) ---
    result = struct( ...
        'config',     struct('NumDel',          numDel, ...
                             'MaxCandidates',   maxSolutions, ...
                             'Values',          values, ...
                             'Sense',           sense, ...
                             'VerbFlag',        verbFlag, ...
                             'SelectedRxnList', {selectedRxnList}), ...
        'biomassRxn', biomassRxn, ...
        'targetRxn',  targetRxn, ...
        'outputFile', outputFile, ...
        'matFile',    '', ...
        'rows',       struct('solutionID', {}, 'setSize',    {}, ...
                             'rxnID',      {}, 'rxnName',    {}, ...
                             'subsystems', {}, 'grRule',     {}), ...
        'targets',    struct('rxnList',         [], ...
                             'appearanceCount', [], ...
                             'wtBiomass',       [], ...
                             'targetMaxFlux',   []));

    % --- Open CSV (file or stdout) and write header ---
    % (No Direction column per MDP-extension)
    headerLine = 'KO_ID,KO_Name,Subsystems,GrRule';
    if output
        outputFile = char(outputFile);
        fid = fopen(outputFile, 'w');
        fprintf(fid, '%s\n', headerLine);
    else
        fprintf('%s\n', headerLine);
    end

    % --- Core loop (logic preserved from original) ---
    options   = struct('targetRxn', targetRxn, 'numDel', numDel);
    constrOpt = struct('rxnList', {{biomassRxn}}, 'values', values*sol.f, 'sense', sense);
    prevSolutions = {};
    noSolution = true;

    for i = 1:maxSolutions
        [optKnockSol, ~] = OptKnock(model, selectedRxnList, options, constrOpt, prevSolutions, verbFlag, solutionFileNameTmp);

        if isempty(optKnockSol) || ~isfield(optKnockSol, 'rxnList') || isempty(optKnockSol.rxnList)
            break
        end

        noSolution = false;
        setSize    = numel(optKnockSol.rxnList);
        rxnList_i  = optKnockSol.rxnList;

        % Map numeric indices to IDs (if needed); otherwise pass through
        if isnumeric(rxnList_i)
            rxnIDs = model.rxns(rxnList_i);
        else
            rxnIDs = rxnList_i;
        end

        % --- Emit each KO as a row, grouped by solution ---
        for k = 1:setSize
            row = struct( ...
                'solutionID', i, ...
                'setSize',    setSize, ...
                'rxnID',      rxnIDs{k}, ...
                'rxnName',    '', ...
                'subsystems', '', ...
                'grRule',     '');

            pos = find(strcmp(model.rxns, rxnIDs{k}), 1);
            if ~isempty(pos)
                if isfield(model, 'rxnNames')
                    row.rxnName = char(model.rxnNames{pos});
                end
                if isfield(model, 'subSystems')
                    ss = model.subSystems{pos};
                    if ~iscell(ss)
                        ss = {ss};
                    end
                    row.subsystems = strjoin(ss, ';');
                end
                if isfield(model, 'grRules')
                    row.grRule = char(model.grRules{pos});
                end
            end

            result.rows(end+1) = row;

            % Write to CSV / stdout (no Direction column)
            if output
                fprintf(fid, '%s,%s,%s,%s\n', ...
                    csvEscape(rxnIDs{k}), csvEscape(row.rxnName), ...
                    csvEscape(row.subsystems), csvEscape(row.grRule));
            else
                fprintf('%s,%s,%s,%s\n', ...
                    csvEscape(rxnIDs{k}), csvEscape(row.rxnName), ...
                    csvEscape(row.subsystems), csvEscape(row.grRule));
            end
        end

        % --- Visual separator between solutions (blank line) ---
        if output
            fprintf(fid, '\n');
        else
            fprintf('\n');
        end

        prevSolutions{end+1} = optKnockSol.rxnList;
    end

    if output
        fclose(fid);
    end

    % --- No-solution fallback: rewrite file with just "# No solutions found" ---
    if noSolution
        if output
            fid = fopen(outputFile, 'w');
            fprintf(fid, '# No solutions found\n');
            fclose(fid);
        else
            fprintf('# No solutions found\n');
        end
    end

    % --- Populate result.targets ---
    result.targets.wtBiomass     = wtBiomass;
    result.targets.targetMaxFlux = targetMax;

    % Flatten all solutions into unique IDs (for KO target summary)
    flatIds = {};
    for s = 1:numel(prevSolutions)
        r = prevSolutions{s};
        if isnumeric(r)
            flatIds = [flatIds; model.rxns(r(:))];
        elseif iscell(r)
            flatIds = [flatIds; r(:)];
        end
    end
    if ~isempty(flatIds)
        [uniqueIds, ~, idxMap] = unique(flatIds);
        counts = accumarray(idxMap, 1);
        result.targets.rxnList         = uniqueIds;
        result.targets.appearanceCount = counts;
    end

    % --- Persist result.mat (mirror algFseof) ---
    if ~isempty(algoDir) && ~exist(algoDir, 'dir')
        mkdir(algoDir);
    end
    result.matFile = matFile;
    save(result.matFile, 'result');
    fprintf('[optknock] Saved `result` to %s\n', result.matFile);
end


function safe = safeTargetName(targetRxn)
    safe = regexprep(char(targetRxn), '[^A-Za-z0-9_-]+', '_');
    safe = regexprep(safe, '^_+|_+$', '');
    if isempty(safe), safe = 'target'; end
end


function s = csvEscape(s)
% csvEscape  Quote a string for CSV output if it contains a comma, double
%   quote, or newline. Embedded double quotes are doubled per RFC 4180.
    if isempty(s), return; end
    if any(s == ',') || any(s == '"') || any(s == sprintf('\n')) || any(s == sprintf('\r'))
        s = strrep(s, '"', '""');
        s = ['"' s '"'];
    end
end
