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
%       .metabolites   - Cell array of metabolite IDs for each reaction
%       .coefficients  - Matrix of stoichiometric coefficients
%
% Data Format (New TSV):
%   Type       RxnName    ModelMetabolite    Coefficient    Cond1    Cond2    ...
%   constraint biomass    {biomass}          {-1}           0.08     0.09    ...
%   constraint EX_glc__D_e {EX_glc__D_e}      {-1}           -10      -12     ...
%   flux       HEX1       {glc__D_c;atp_c;g6p_c} {1;1;-1}   0.91     1.19    ...
%   flux       PGI        {g6p_c;f6p_c}      {1;-1}        0.65     0.85    ...
%
%   - constraint rows: define model bounds and objectives
%   - flux rows: define central carbon metabolism fluxes with stoichiometry

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

    % Find condition columns (Cond1, Cond2, ...)
    varNames = data.Properties.VariableNames;
    condCols = find(cellfun(@(x) startsWith(x, 'Cond'), varNames));
    conditionNames = varNames(condCols);
    nCond = length(condCols);

    % Initialize output structure
    C13Data.constraints = struct();
    C13Data.reactions = {};
    C13Data.conditions = conditionNames;
    C13Data.fluxes = [];
    C13Data.metabolites = {};
    C13Data.coefficients = [];

    % Get data as arrays
    types = data{:, 'Type'};
    rxnNames = data{:, 'RxnName'};
    modelMets = data{:, 'ModelMetabolite'};
    coeffs = data{:, 'Coefficient'};

    % Separate constraint and flux rows
    constraintIdx = find(strcmp(types, 'constraint'));
    fluxIdx = find(strcmp(types, 'flux'));

    fprintf('[load13CData] Found %d constraints and %d flux reactions\n', ...
        length(constraintIdx), length(fluxIdx));

    % --- Parse constraint rows ---
    for i = 1:length(constraintIdx)
        idx = constraintIdx(i);
        rxn = rxnNames{idx};

        % Parse ModelMetabolite (e.g., '{biomass}' or '{EX_glc__D_e}')
        metStr = modelMets{idx};
        if ischar(metStr)
            metStr = strtrim(metStr);
            % Remove curly braces
            metStr = strrep(metStr, '{', '');
            metStr = strrep(metStr, '}', '');
        end

        % Parse Coefficient (e.g., '{-1}')
        coefStr = coeffs{idx};
        if ischar(coefStr)
            coefStr = strtrim(coefStr);
            coefStr = strrep(coefStr, '{', '');
            coefStr = strrep(coefStr, '}', '');
            coefVal = str2double(strsplit(coefStr, ';'));
        else
            coefVal = coefStr;
        end

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
    end

    % --- Parse flux rows ---
    fluxData = [];
    metabolitesData = {};
    coefficientsData = [];

    for i = 1:length(fluxIdx)
        idx = fluxIdx(i);
        rxn = rxnNames{idx};

        % Parse ModelMetabolite (e.g., '{glc__D_c;atp_c;g6p_c}')
        metStr = modelMets{idx};
        if ischar(metStr)
            metStr = strtrim(metStr);
            metStr = strrep(metStr, '{', '');
            metStr = strrep(metStr, '}', '');
            mets = strsplit(metStr, ';');
            mets = mets(~cellfun(@isempty, mets));
        else
            mets = {metStr};
        end

        % Parse Coefficient (e.g., '{1;1;-1}')
        coefStr = coeffs{idx};
        if ischar(coefStr)
            coefStr = strtrim(coefStr);
            coefStr = strrep(coefStr, '{', '');
            coefStr = strrep(coefStr, '}', '');
            coefVals = str2double(strsplit(coefStr, ';'));
            coefVals = coefVals(~isnan(coefVals));
        else
            coefVals = coefStr;
        end

        % Get flux values for each condition
        fluxVals = NaN(1, nCond);
        for c = 1:nCond
            val = data{idx, condCols(c)};
            if ~ismissing(val) && ~isnan(val)
                fluxVals(c) = val;
            end
        end

        % Store
        C13Data.reactions{end+1} = rxn;
        metabolitesData{end+1} = mets;
        coefficientsData{end+1, 1} = coefVals;
        fluxData = [fluxData; fluxVals]; %#ok<AGROW>
    end

    C13Data.fluxes = fluxData;
    C13Data.metabolites = metabolitesData;
    C13Data.coefficients = coefficientsData;

    % Print sample flux reactions
    fprintf('[load13CData] Sample flux reactions:\n');
    for i = 1:min(5, length(C13Data.reactions))
        mets = C13Data.metabolites{i};
        coefs = C13Data.coefficients{i};
        fprintf('  %s: %s (%s)\n', C13Data.reactions{i}, ...
            strjoin(mets(1:min(3, end)), ','), ...
            mat2str(coefs(1:min(3, end)), 2));
    end

    % Match reactions to model if provided
    if ~isempty(model)
        C13Data = matchReactionsToModel(C13Data, model);
    end
end

%% ------------------- Model Matching -------------------

function C13Data = matchReactionsToModel(C13Data, model)
% Match experimental reactions to model reactions

    if isfield(model, 'rxns')
        modelRxns = model.rxns;
    else
        warning('Model does not have rxns field');
        return;
    end

    % Match each reaction
    matchedIdx = zeros(length(C13Data.reactions), 1);
    for i = 1:length(C13Data.reactions)
        expRxn = C13Data.reactions{i};
        idx = find(strcmp(modelRxns, expRxn));
        if ~isempty(idx)
            matchedIdx(i) = idx(1);
        end
    end

    nMatched = sum(matchedIdx > 0);
    fprintf('[load13CData] Matched %d/%d reactions to model\n', nMatched, length(C13Data.reactions));

    % Store matching info
    C13Data.matchedIdx = matchedIdx;
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
end
