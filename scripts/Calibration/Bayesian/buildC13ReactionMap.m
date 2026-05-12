function [C13Map, validationReport] = buildC13ReactionMap(C13Data, ecModel)
% buildC13ReactionMap
%   Builds a mapping between 13C flux data reactions and ecModel reactions.
%
%   Logic:
%   1. User's reaction name (e.g., 'r_0534' or '{r_0658;r_2131}') is cleaned
%   2. Cleaned name is matched to model reactions (via cleanedToIndices)
%   3. All matching model reactions are collected (for split/isozyme reactions)
%   4. Final direction = user_direction × model_reaction_direction
%
% Input:
%   C13Data  - 13C flux data structure from load13CData
%              Required fields: reactions, types, directions, nCarbon
%   ecModel  - Enzyme-constrained model
%
% Output:
%   C13Map - Structure containing reaction mappings:
%       .reactions       - User's reaction names
%       .modelIndices    - Model reaction indices (may be multiple for split reactions)
%       .directions     - Combined directions (user × model)
%       .nCarbon        - Carbon count
%       .matchStatus    - 'matched', 'unmatched', 'multi_matched'
%
%   validationReport - Validation statistics

    %% ---- Step 1: Initialize and detect model type ----
    fprintf('[buildC13ReactionMap] Starting reaction mapping construction...\n');

    % Get model type
    ecModelType = [];
    if isfield(ecModel, 'enzymeConstraints') && isfield(ecModel.enzymeConstraints, 'ecModeltype')
        ecModelType = ecModel.enzymeConstraints.ecModeltype;
    end
    isBasicModel = strcmp(ecModelType, 'basic');

    fprintf('[buildC13ReactionMap] Model type: %s\n', mat2str(ecModelType));

    % Get model reaction info
    modelRxns = ecModel.rxns;
    nModelRxns = length(modelRxns);

    % Build cleaned model reaction names and store suffix info
    cleanedModelRxns = cell(nModelRxns, 1);
    suffixInfoList = cell(nModelRxns, 1);  % Store suffix info for direction detection

    for i = 1:nModelRxns
        [cleanName, suffixInfo] = cleanReactionName(modelRxns{i}, ecModelType);
        cleanedModelRxns{i} = cleanName;
        suffixInfoList{i} = suffixInfo;
    end

    % Build reverse index: cleaned name -> ALL original indices
    uniqueCleanedNames = unique(cleanedModelRxns);
    cleanedToIndices = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:length(uniqueCleanedNames)
        cleanedToIndices(uniqueCleanedNames{i}) = [];
    end
    for i = 1:nModelRxns
        cleanedToIndices(cleanedModelRxns{i}) = [cleanedToIndices(cleanedModelRxns{i}), i];
    end

    %% ---- Step 2: Get 13C data ----
    c13Reactions = C13Data.reactions;  % Cell array: single string or {r1;r2;r3}
    c13Directions = C13Data.directions; % User-provided directions
    c13Types = C13Data.types;          % 'flux' or 'constraint'
    c13Carbon = C13Data.nCarbon;        % Carbon count
    nC13Reactions = length(c13Reactions);

    fprintf('[buildC13ReactionMap] Processing %d 13C reactions\n', nC13Reactions);

    % Count types
    nConstraint = sum(strcmp(c13Types, 'constraint'));
    nFlux = sum(strcmp(c13Types, 'flux'));
    fprintf('[buildC13ReactionMap]   - %d constraint reactions\n', nConstraint);
    fprintf('[buildC13ReactionMap]   - %d flux reactions\n', nFlux);

    %% ---- Step 3: Build reaction mapping ----
    C13Map = struct();
    C13Map.reactions = cell(nC13Reactions, 1);
    C13Map.originalNames = c13Reactions;
    C13Map.types = c13Types;
    C13Map.modelIndices = cell(nC13Reactions, 1);  % May contain multiple indices for split reactions
    C13Map.directions = cell(nC13Reactions, 1);   % Combined directions
    C13Map.nCarbon = zeros(nC13Reactions, 1);
    C13Map.matchStatus = cell(nC13Reactions, 1);
    C13Map.warnings = {};

    % Validation tracking
    validationReport = struct();
    validationReport.totalReactions = nC13Reactions;
    validationReport.matchedCount = 0;
    validationReport.unmatchedCount = 0;
    validationReport.multiMatchedCount = 0;
    validationReport.warnings = {};

    unmatchedList = {};

    %% ---- Step 4: Match each 13C reaction to model reactions ----
    for i = 1:nC13Reactions
        c13Rxn = c13Reactions{i};
        c13Dir = c13Directions{i};  % User-provided direction array
        rxnType = c13Types{i};

        % Parse user's reaction input
        if iscell(c13Rxn)
            % Grouped reaction: {r1;r2;r3}
            userRxnList = c13Rxn;
        else
            % Single reaction
            userRxnList = {c13Rxn};
        end

        % Collect all matched model indices and their inherent directions
        allMatchedIndices = [];
        allModelDirs = [];  % Model's inherent direction (1=forward, 0=reverse)

        for j = 1:length(userRxnList)
            userRxnName = userRxnList{j};

            % Clean the reaction name
            [cleanC13Name, ~] = cleanReactionName(userRxnName, ecModelType);

            % Find ALL matching model reactions
            if isKey(cleanedToIndices, cleanC13Name)
                matchedModelIdx = cleanedToIndices(cleanC13Name);

                % Get model direction for each matched reaction
                for k = 1:length(matchedModelIdx)
                    modelIdx = matchedModelIdx(k);
                    suffixInfo = suffixInfoList{modelIdx};

                    % Get model reaction's inherent direction
                    if isBasicModel
                        % Basic model: prefix determines direction
                        modelDir = (suffixInfo.copyNumber == 1);
                    else
                        % Isozyme/integrated: _REV determines direction
                        modelDir = ~suffixInfo.isReverse;  % true=forward=1, false=reverse=0
                    end

                    allMatchedIndices = [allMatchedIndices, modelIdx];
                    allModelDirs = [allModelDirs, modelDir];
                end
            end
        end

        nMatched = length(allMatchedIndices);

        % Get user's direction for each reaction in the group
        % User direction: 1 = forward, -1 = reverse
        % If user provides fewer directions than matched reactions, default to 1
        userDirs = ones(1, nMatched);
        if ~isempty(c13Dir)
            for j = 1:min(length(c13Dir), nMatched)
                userDirs(j) = c13Dir(j);
            end
        end

        % Combine user direction with model direction
        % - user=1 (forward): EXP contribute +flux, REV contribute +|flux|
        % - user=-1 (reverse): EXP contribute -flux, REV contribute -|flux|
        %
        % In FBA:
        % - EXP reactions have positive flux when running forward
        % - REV reactions have negative flux when running forward (in their direction)
        %
        % Contribution = user_dir × flux for EXP
        % Contribution = -user_dir × flux for REV
        % (because REV flux is negative when running in REV direction)
        %
        % So final direction factor:
        % - EXP: user_dir × 1 = user_dir
        % - REV: user_dir × -1 = -user_dir
        finalDirs = zeros(1, nMatched);
        for k = 1:nMatched
            if allModelDirs(k) == 1
                % Forward in model: use user direction
                finalDirs(k) = userDirs(k);
            else
                % Reverse in model: negate user direction
                finalDirs(k) = -userDirs(k);
            end
        end

        % Store in C13Map
        C13Map.reactions{i} = c13Rxn;
        C13Map.modelIndices{i} = allMatchedIndices;
        C13Map.directions{i} = finalDirs;

        % For flux: use user-provided carbon count
        % For constraint: get from model (ecModel.excarbon)
        if strcmp(rxnType, 'constraint')
            % Get carbon from ecModel.excarbon using first matched index
            if ~isempty(allMatchedIndices) && isfield(ecModel, 'excarbon')
                firstIdx = allMatchedIndices(1);
                if firstIdx <= length(ecModel.excarbon)
                    C13Map.nCarbon(i) = abs(ecModel.excarbon(firstIdx));
                else
                    C13Map.nCarbon(i) = 0;
                end
            else
                C13Map.nCarbon(i) = 0;
            end
        else
            C13Map.nCarbon(i) = c13Carbon(i);  % Use user-provided carbon for flux
        end

        % Determine match status
        if nMatched == 0
            C13Map.matchStatus{i} = 'unmatched';
            validationReport.unmatchedCount = validationReport.unmatchedCount + 1;
            if iscell(c13Rxn)
                unmatchedList{end+1, 1} = ['{', strjoin(c13Rxn, ';'), '}']; %#ok<AGROW>
            else
                unmatchedList{end+1, 1} = c13Rxn; %#ok<AGROW>
            end
        elseif nMatched == 1
            C13Map.matchStatus{i} = 'matched';
            validationReport.matchedCount = validationReport.matchedCount + 1;
        else
            C13Map.matchStatus{i} = 'multi_matched';
            validationReport.multiMatchedCount = validationReport.multiMatchedCount + 1;
            validationReport.matchedCount = validationReport.matchedCount + 1;
        end
    end

    %% ---- Step 5: Print summary ----
    fprintf('\n[buildC13ReactionMap] Mapping Summary:\n');
    fprintf('  Total 13C reactions: %d\n', validationReport.totalReactions);
    fprintf('    - Constraint reactions: %d\n', nConstraint);
    fprintf('    - Flux reactions: %d\n', nFlux);
    fprintf('  Matched (single):   %d\n', validationReport.matchedCount - validationReport.multiMatchedCount);
    fprintf('  Multi-matched:      %d\n', validationReport.multiMatchedCount);
    fprintf('  Unmatched:          %d\n', validationReport.unmatchedCount);

    if ~isempty(unmatchedList)
        fprintf('\n[buildC13ReactionMap] Unmatched reactions:\n');
        for i = 1:min(10, length(unmatchedList))
            fprintf('  - %s\n', unmatchedList{i});
        end
        if length(unmatchedList) > 10
            fprintf('  ... and %d more\n', length(unmatchedList) - 10);
        end
    end

    % Print multi-matched reactions
    multiMatchedReactions = {};
    for i = 1:nC13Reactions
        if strcmp(C13Map.matchStatus{i}, 'multi_matched')
            rxn = C13Map.reactions{i};
            if iscell(rxn)
                multiMatchedReactions{end+1, 1} = ['{', strjoin(rxn, ';'), '}']; %#ok<AGROW>
            else
                multiMatchedReactions{end+1, 1} = rxn; %#ok<AGROW>
            end
        end
    end

    if ~isempty(multiMatchedReactions)
        fprintf('\n[buildC13ReactionMap] Multi-matched reactions (split/isozyme):\n');
        for i = 1:min(10, length(multiMatchedReactions))
            fprintf('  - %s\n', multiMatchedReactions{i});
        end
        if length(multiMatchedReactions) > 10
            fprintf('  ... and %d more\n', length(multiMatchedReactions) - 10);
        end
    end

    validationReport.warnings = C13Map.warnings;
end

%% ============== Helper Functions ==============

function [cleanName, suffixInfo] = cleanReactionName(rxnName, ecModelType)
% cleanReactionName
%   Removes model-specific suffixes/prefixes from reaction names
%
%   Naming conventions:
%   - basic model: 001_RXN, 002_RXN, ... (prefix)
%   - integrated/isozyme:
%       - RXN_EXP_1, RXN_EXP_2 (isozyme only, forward)
%       - RXN_REV_EXP_1, RXN_REV_EXP_2 (isozyme, reverse)
%       - RXN_f, RXN_b (reversible split)
%       - RXN_f_EXP_1, RXN_b_EXP_1 (combined)
%       - RXN_forward, RXN_reverse (COBRA style)

    suffixInfo = struct();
    suffixInfo.pattern = '';
    suffixInfo.type = '';
    suffixInfo.isReverse = false;
    suffixInfo.copyNumber = NaN;

    cleanName = rxnName;

    % === Check for prefix patterns first (basic model: 001_, 002_, ...) ===
    prefixMatch = regexp(cleanName, '^(\d{3,})_(.+)$', 'tokens');
    if ~isempty(prefixMatch)
        cleanName = prefixMatch{1}{2};
        suffixInfo.type = 'prefix';
        suffixInfo.pattern = [prefixMatch{1}{1}, '_'];
        suffixInfo.copyNumber = str2double(prefixMatch{1}{1});
        return;
    end

    % === Loop through suffix patterns ===
    while true
        removed = false;

        % _REV (reversible marker)
        if endsWith(cleanName, '_REV')
            cleanName = cleanName(1:end-4);
            suffixInfo.type = 'suffix';
            suffixInfo.pattern = [suffixInfo.pattern, '_REV'];
            suffixInfo.isReverse = true;
            removed = true;
            continue;
        end

        % _EXP_1, _EXP_2, ... (isozyme expansion)
        expMatch = regexp(cleanName, '_EXP_(\d+)$', 'tokens');
        if ~isempty(expMatch)
            cleanName = regexprep(cleanName, '_EXP_\d+$', '');
            suffixInfo.type = 'suffix';
            suffixInfo.pattern = [suffixInfo.pattern, '_EXP_', expMatch{1}{1}];
            removed = true;
            continue;
        end

        % _forward, _reverse (COBRA style)
        if endsWith(cleanName, '_forward')
            cleanName = cleanName(1:end-8);
            suffixInfo.type = 'suffix';
            suffixInfo.pattern = [suffixInfo.pattern, '_forward'];
            removed = true;
            continue;
        end
        if endsWith(cleanName, '_reverse')
            cleanName = cleanName(1:end-8);
            suffixInfo.type = 'suffix';
            suffixInfo.pattern = [suffixInfo.pattern, '_reverse'];
            suffixInfo.isReverse = true;
            removed = true;
            continue;
        end

        % _f, _b (irreversible split)
        if endsWith(cleanName, '_f') && ~endsWith(cleanName, '_cf') && ~endsWith(cleanName, '_uf')
            cleanName = cleanName(1:end-2);
            suffixInfo.type = 'suffix';
            suffixInfo.pattern = [suffixInfo.pattern, '_f'];
            removed = true;
            continue;
        end
        if endsWith(cleanName, '_b') && ~endsWith(cleanName, '_cb') && ~endsWith(cleanName, '_ub')
            cleanName = cleanName(1:end-2);
            suffixInfo.type = 'suffix';
            suffixInfo.pattern = [suffixInfo.pattern, '_b'];
            suffixInfo.isReverse = true;
            removed = true;
            continue;
        end

        % No more suffixes to remove
        break;
    end
end
