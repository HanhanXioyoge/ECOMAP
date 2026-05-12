function C13Data = load13CData(filePath, model)
% load13CData
%   Reads 13C flux data from a TSV file and prepares it for RMSE calculation.
%
% Input:
%   filePath    - Path to the 13C flux data TSV file
%   model       - ecModel (used for reaction matching, optional)
%
% Output:
%   C13Data    - Structure containing:
%       .constraints  - Struct with model constraints (lb, ub, objective)
%       .reactions     - Cell array of reaction names (flux rows only)
%       .conditions    - Cell array of condition names (Cond1, Cond2, ...)
%       .fluxes        - Matrix of flux values (nReactions x nConditions)
%       .directions    - Cell array of direction arrays for each flux reaction
%       .nCarbon       - User-provided carbon count for flux reactions
%
% Data Format (Simplified TSV):
%   Type       RxnName              Carbon    Direction    Cond1    Cond2    ...
%   constraint r_2111                                    0.405    0.150    ...
%   constraint r_1714                                   -16.731   -1.560   ...
%   flux       r_0534              6        {1}        16.731    1.560    ...
%   flux       {r_0658;r_2131}     6        {1;1}       1.008    1.140    ...
%   flux       {r_0454;r_1021;r_0452}  4     {-1;1;-1}   0.102    0.989    ...
%
%   - constraint rows: define model bounds and objectives (Direction is ignored)
%   - flux rows: define central carbon metabolism fluxes
%     - RxnName: single reaction ID or {r1;r2;r3} for multiple reactions
%     - Direction: {1} for single, {1;-1;1} for multiple (1=forward, -1=reverse)
%     - Carbon: carbon atoms in the substrate/product

    if nargin < 2
        model = [];
    end

    % Read the file
    if ~exist(filePath, 'file')
        error('13C data file not found: %s', filePath);
    end

    fprintf('[load13CData] Loading 13C data from: %s\n', filePath);

    % Read TSV file
    rawData = readtable(filePath, 'FileType', 'text', 'ReadVariableNames', true);

    % Get column names
    varNames = rawData.Properties.VariableNames;

    % Check for new format (has Type column)
    if ismember('Type', varNames)
        C13Data = parseNewFormat(rawData, model);
    else
        % Fallback to old format detection
        C13Data = parseOldFormat(rawData);
    end

    fprintf('[load13CData] Loaded %d flux reactions x %d conditions\n', ...
        length(C13Data.reactions), length(C13Data.conditions));
    fprintf('[load13CData] Loaded %d constraints\n', length(fieldnames(C13Data.constraints)));
end

%% ------------------- New Format Parser -------------------

function C13Data = parseNewFormat(data, model)
% Parse new format with Type column
%
% Expected columns: Type, RxnName, Carbon, Direction, Cond1, Cond2, ...

    % Find condition columns (Cond1, Cond2, ...)
    varNames = data.Properties.VariableNames;
    condCols = find(cellfun(@(x) startsWith(x, 'Cond'), varNames));
    conditionNames = varNames(condCols);
    nCond = length(condCols);

    % Initialize output structure
    C13Data.constraints = struct();
    C13Data.reactions = {};      % Cell array of reaction names (single or cell array for grouped)
    C13Data.types = {};          % Store reaction types ('flux' or 'constraint')
    C13Data.conditions = conditionNames;
    C13Data.fluxes = [];
    C13Data.directions = {};     % Direction arrays for flux reactions
    C13Data.nCarbon = [];       % User-provided carbon count for flux reactions

    % Get data as arrays
    types = data{:, 'Type'};
    rxnNames = data{:, 'RxnName'};
    carbonVals = data{:, 'Carbon'};

    % Check if Direction column exists
    hasDirectionCol = ismember('Direction', varNames);
    if hasDirectionCol
        directionVals = data{:, 'Direction'};
    else
        directionVals = cell(height(data), 1);
    end

    % Separate constraint and flux rows
    constraintIdx = find(strcmp(types, 'constraint'));
    fluxIdx = find(strcmp(types, 'flux'));

    fprintf('[load13CData] Found %d constraints and %d flux reactions\n', ...
        length(constraintIdx), length(fluxIdx));

    % --- Parse constraint rows ---
    for i = 1:length(constraintIdx)
        idx = constraintIdx(i);
        rxn = rxnNames{idx};

        % Get flux values for each condition
        fluxVals = NaN(1, nCond);
        for c = 1:nCond
            val = data{idx, condCols(c)};
            if ~ismissing(val) && ~isnan(val)
                fluxVals(c) = val;
            end
        end

        % Determine constraint type based on reaction name
        rxnLower = lower(rxn);

        % Store constraint based on reaction type
        if contains(rxnLower, 'biomass')
            % Biomass objective - store as objective
            C13Data.constraints.biomass = struct('objective', fluxVals);
            fprintf('  Constraint: %s (objective), values: %s\n', rxn, mat2str(fluxVals(~isnan(fluxVals)), 4));
        elseif startsWith(rxn, 'EX_')
            % Exchange reaction - store as bounds
            C13Data.constraints.(rxn) = struct('lb', fluxVals, 'ub', zeros(1, nCond) + 1000);
            fprintf('  Constraint: %s (exchange), lb: %s\n', rxn, mat2str(fluxVals(~isnan(fluxVals)), 4));
        else
            % Other constraints
            C13Data.constraints.(rxn) = struct('lb', fluxVals, 'ub', zeros(1, nCond) + 1000);
            fprintf('  Constraint: %s, lb: %s\n', rxn, mat2str(fluxVals(~isnan(fluxVals)), 4));
        end

        % Add to main data arrays (for RMSE calculation)
        C13Data.reactions{end+1} = rxn;
        C13Data.types{end+1} = 'constraint';
        C13Data.directions{end+1} = [];  % No direction for constraint
        C13Data.nCarbon(end+1, 1) = NaN;  % Constraint gets carbon from ecModel.excarbon
        C13Data.fluxes = [C13Data.fluxes; fluxVals]; %#ok<AGROW>
    end

    % --- Parse flux rows ---
    for i = 1:length(fluxIdx)
        idx = fluxIdx(i);
        rxnRaw = rxnNames{idx};

        % Parse RxnName: single reaction or {r1;r2;r3}
        [rxnList, isMulti] = parseReactionName(rxnRaw);

        % Parse Direction: {1} or {1;-1;1}
        dirStr = '';
        if hasDirectionCol
            dirStr = directionVals{idx};
        end
        dirArray = parseDirections(dirStr, length(rxnList));

        % Get user-provided carbon count
        if iscell(carbonVals)
            cVal = carbonVals{idx};
            if ischar(cVal)
                cVal = str2double(strrep(strrep(cVal, '{', ''), '}', ''));
            end
        else
            cVal = carbonVals(idx);
        end

        % Get flux values for each condition
        fluxVals = NaN(1, nCond);
        for c = 1:nCond
            val = data{idx, condCols(c)};
            if ~ismissing(val) && ~isnan(val)
                fluxVals(c) = val;
            end
        end

        % Store in main data arrays
        C13Data.reactions{end+1} = rxnList;  % Cell array of reaction names
        C13Data.types{end+1} = 'flux';
        C13Data.directions{end+1} = dirArray;  % Direction array
        C13Data.nCarbon(end+1, 1) = cVal;  % User-provided carbon count
        C13Data.fluxes = [C13Data.fluxes; fluxVals]; %#ok<AGROW>

        % Print info for multi-reaction fluxes
        if isMulti
            dirParts = arrayfun(@(x) sprintf('%.0f', x), dirArray, 'UniformOutput', false);
            fprintf('  Flux: {%s} x %d carbons, directions: {%s}\n', ...
                strjoin(rxnList, ';'), cVal, strjoin(dirParts, ';'));
        end
    end

    % Print sample flux reactions
    fprintf('[load13CData] Sample flux reactions:\n');
    count = 0;
    for i = 1:length(C13Data.reactions)
        if strcmp(C13Data.types{i}, 'flux')
            rxn = C13Data.reactions{i};
            nC = C13Data.nCarbon(i);
            dirs = C13Data.directions{i};
            if iscell(rxn)
                rxnStr = ['{', strjoin(rxn, ';'), '}'];
            else
                rxnStr = rxn;
            end
            % Format directions properly
            if isempty(dirs)
                dirStr = '{}';
            else
                dirParts = arrayfun(@(x) sprintf('%.0f', x), dirs, 'UniformOutput', false);
                dirStr = ['{', strjoin(dirParts, ';'), '}'];
            end
            fprintf('  %s: C=%d, dir=%s\n', rxnStr, nC, dirStr);
            count = count + 1;
            if count >= 5
                break;
            end
        end
    end

    % Match reactions to model if provided
    if ~isempty(model)
        C13Data = matchReactionsToModel(C13Data, model);
    end
end

%% ------------------- Helper Functions -------------------

function [rxnList, isMulti] = parseReactionName(rxnRaw)
% Parse reaction name that can be:
%   - Single: 'r_0534'
%   - Multiple: '{r_0658;r_2131}'

    rxnRaw = strtrim(rxnRaw);
    if startsWith(rxnRaw, '{')
        % Multi-reaction format: {r1;r2;r3}
        rxnRaw = strrep(rxnRaw, '{', '');
        rxnRaw = strrep(rxnRaw, '}', '');
        rxnList = strsplit(rxnRaw, ';');
        rxnList = rxnList(~cellfun(@isempty, rxnList));
        rxnList = strtrim(rxnList);
        isMulti = true;
    else
        % Single reaction
        rxnList = {rxnRaw};
        isMulti = false;
    end
end

function dirArray = parseDirections(dirStr, nExpected)
% Parse direction string that can be:
%   - Empty: '' (for single reaction, default to 1)
%   - Single: '{1}' or '{-1}' (expand to nExpected if multi)
%   - Multiple: '{1;-1;1}'

    dirArray = ones(1, nExpected);  % Default to forward

    % Handle cell array input
    if iscell(dirStr)
        if isempty(dirStr)
            return;
        end
        dirStr = dirStr{1};
        if iscell(dirStr)
            return;  % Still cell, use default
        end
    end

    % Handle empty char
    if ~ischar(dirStr)
        return;  % Not a string, use default
    end

    dirStr = strtrim(dirStr);
    if isempty(dirStr)
        return;
    end
    if startsWith(dirStr, '{')
        dirStr = strrep(dirStr, '{', '');
        dirStr = strrep(dirStr, '}', '');
        parts = strsplit(dirStr, ';');
        parts = parts(~cellfun(@isempty, parts));
        dirVals = str2double(parts);

        if length(dirVals) == 1 && nExpected > 1
            % Single direction for multi-reaction, apply to all
            dirArray(:) = dirVals;
        else
            % Match length
            nUse = min(length(dirVals), nExpected);
            dirArray(1:nUse) = dirVals(1:nUse);
        end
    else
        % Plain number
        val = str2double(dirStr);
        if ~isnan(val)
            dirArray(:) = val;
        end
    end
end

%% ------------------- Model Matching -------------------

function C13Data = matchReactionsToModel(C13Data, model)
% Match experimental reactions to model reactions
% For multi-reaction fluxes, this just stores the names

    if isfield(model, 'rxns')
        modelRxns = model.rxns;
    else
        warning('Model does not have rxns field');
        return;
    end

    % For each reaction, store the reaction list
    % (Actual matching is done in buildC13ReactionMap)
    nRxns = length(C13Data.reactions);
    matchedCount = 0;
    unmatchedReactions = {};

    for i = 1:nRxns
        rxn = C13Data.reactions{i};

        if strcmp(C13Data.types{i}, 'flux')
            % For flux, check if reactions exist in model
            if iscell(rxn)
                % Multi-reaction
                allFound = true;
                for j = 1:length(rxn)
                    if isempty(find(strcmp(modelRxns, rxn{j}), 1))
                        allFound = false;
                        break;
                    end
                end
                if allFound
                    matchedCount = matchedCount + 1;
                else
                    unmatchedReactions{end+1} = ['{', strjoin(rxn, ';'), '}']; %#ok<AGROW>
                end
            else
                % Single reaction
                if ~isempty(find(strcmp(modelRxns, rxn), 1))
                    matchedCount = matchedCount + 1;
                else
                    unmatchedReactions{end+1} = rxn; %#ok<AGROW>
                end
            end
        end
    end

    fprintf('[load13CData] Matched %d/%d flux reactions to model\n', matchedCount, sum(strcmp(C13Data.types, 'flux')));
    if ~isempty(unmatchedReactions)
        fprintf('[load13CData] Unmatched reactions (first 10):\n');
        for i = 1:min(10, length(unmatchedReactions))
            fprintf('  - %s\n', unmatchedReactions{i});
        end
    end

    C13Data.modelRxns = modelRxns;
end

%% ------------------- Old Format Parser (Legacy) -------------------

function C13Data = parseOldFormat(data)
% Parse old format (for backward compatibility)

    varNames = data.Properties.VariableNames;
    nVars = length(varNames);

    if nVars == 3
        % Long format
        conditions = unique(data{:, 1});
        reactions = unique(data{:, 2});
        nCond = length(conditions);
        nRxn = length(reactions);

        fluxes = NaN(nRxn, nCond);
        for i = 1:height(data)
            condIdx = find(strcmp(conditions, data{i, 1}));
            rxnIdx = find(strcmp(reactions, data{i, 2}));
            fluxes(rxnIdx, condIdx) = data{i, 3};
        end

        C13Data.reactions = reactions;
        C13Data.conditions = conditions;
        C13Data.fluxes = fluxes;
        C13Data.format = 'long';
    else
        % Wide format
        reactions = data{:, 1};
        condNames = varNames(2:end);
        fluxes = table2array(data(:, 2:end));
        fluxes = str2double(fluxes);

        C13Data.reactions = reactions;
        C13Data.conditions = condNames;
        C13Data.fluxes = fluxes;
        C13Data.format = 'wide';
    end

    % Empty constraints for old format
    C13Data.constraints = struct();
    C13Data.directions = {};
    C13Data.nCarbon = [];
end
