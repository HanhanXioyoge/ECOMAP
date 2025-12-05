function OUT = merge_match_tables(matchdata, DeepLearningModel, varargin)
% MERGE_MATCH_TABLES
%   Merge tables from multiple saved "match" .mat files by deep-learning model name.
%
%   OUT = merge_match_tables(matchdata, DeepLearningModel)
%   OUT = merge_match_tables(matchdata, DeepLearningModel, 'AddSource', true, 'SourceVarName','source')
%
% INPUT
%   matchdata         : cellstr, e.g. {'eciML1515_match.mat','ecYeast_match.mat', ...}
%                       Each .mat is expected to contain a struct (commonly named 'MATCH')
%                       with fields like 'DLKcat','UniKP','CatPred', and each field is a table.
%   DeepLearningModel : cellstr, subset of model names to merge, e.g. {'DLKcat','UniKP','CatPred'}
%
% NAME-VALUE PAIRS (optional)
%   'AddSource'       : logical (default: false). If true, add a column indicating source file name.
%   'SourceVarName'   : char/string (default: 'source_file'). Name of the source column if AddSource=true.
%
% OUTPUT
%   OUT               : struct. Fields are model names in DeepLearningModel. Each field is a merged table.
%
% BEHAVIOR & ASSUMPTIONS
%   - A .mat file may contain a variable named 'MATCH'. If not, the function searches for a struct variable
%     whose fields include at least one of DeepLearningModel (best-effort).
%   - A model field may be missing in a given file; that file is skipped for that model (warning only).
%   - Tables are vertically concatenated. If their variable sets differ, the function creates missing columns
%     with appropriate "missing" values (NaN, "", NaT, etc.) to make schemas consistent before concatenation.
%   - No de-duplication is performed (you can add it later if needed).
%
% EXAMPLE
%   matchdata = {'eciML1515_match.mat','ecYeast_match.mat','ecHuman_match.mat','eciCW773_match.mat'};
%   models    = {'DLKcat','UniKP','CatPred'};
%   OUT = merge_match_tables(matchdata, models, 'AddSource', true);
%
%   % OUT.DLKcat, OUT.UniKP, OUT.CatPred are merged tables (if requested & available).
%
% Author: (your name), 2025-11-06

    % ----------------------------- Parse inputs -----------------------------
    if nargin < 2
        error('merge_match_tables requires (matchdata, DeepLearningModel).');
    end
    if ~iscellstr(matchdata) && ~all(cellfun(@(x)ischar(x)||isstring(x), matchdata))
        error('matchdata must be a cell array of char/string file names.');
    end
    if ~iscellstr(DeepLearningModel) && ~all(cellfun(@(x)ischar(x)||isstring(x), DeepLearningModel))
        error('DeepLearningModel must be a cell array of char/string model names.');
    end
    matchdata = cellstr(string(matchdata));
    DeepLearningModel = cellstr(string(DeepLearningModel));

    ip = inputParser;
    ip.addParameter('AddSource', false, @(x)islogical(x)&&isscalar(x));
    ip.addParameter('SourceVarName', 'source_file', @(x)ischar(x)||isstring(x));
    ip.parse(varargin{:});
    addSource = ip.Results.AddSource;
    sourceVar = char(ip.Results.SourceVarName);

    % Prepare output as struct of empty tables
    OUT = struct();
    for k = 1:numel(DeepLearningModel)
        OUT.(DeepLearningModel{k}) = table(); %#ok<STRNU>
    end

    % For collecting per-model tables before final concat
    buckets = containers.Map();
    for k = 1:numel(DeepLearningModel)
        buckets(DeepLearningModel{k}) = {}; % cell array of tables
    end

    % ----------------------------- Load & collect ---------------------------
    for i = 1:numel(matchdata)
        f = matchdata{i};
        f = fullfile(findECOMAProot, 'DLmode_evaluation','data', f);
        if ~isfile(f)
            warning('merge_match_tables:MissingFile', 'File not found: %s (skipped)', f);
            continue;
        end

        S = load(f);
        try
            M = extractMatchStruct(S, DeepLearningModel);
        catch ME
            warning('merge_match_tables:BadStruct', ...
                'File %s does not contain a usable MATCH-like struct. (%s) Skipped.', f, ME.message);
            continue;
        end

        for k = 1:numel(DeepLearningModel)
            mdl = DeepLearningModel{k};
            if ~isfield(M, mdl)
                warning('merge_match_tables:MissingModel', ...
                    'File %s: struct has no field "%s" (skipped for this model).', f, mdl);
                continue;
            end
            T = M.(mdl);
            if ~istable(T)
                warning('merge_match_tables:NotATable', ...
                    'File %s: field "%s" is not a table (skipped for this model).', f, mdl);
                continue;
            end
            % Optionally add source column (as string)
            if addSource
                srcCol = repmat(string(f), height(T), 1);
                % If there is already a col named sourceVar, avoid clash by appending
                if any(strcmp(T.Properties.VariableNames, sourceVar))
                    altName = genUniqueVarName(T.Properties.VariableNames, [sourceVar '_2']);
                    T.(altName) = srcCol;
                else
                    T.(sourceVar) = srcCol;
                end
            end
            lst = buckets(mdl);
            lst{end+1} = T; %#ok<AGROW>
            buckets(mdl) = lst;
        end
    end

    % ----------------------------- Merge by model ---------------------------
    for k = 1:numel(DeepLearningModel)
        mdl = DeepLearningModel{k};
        tlist = buckets(mdl);
        if isempty(tlist)
            % leave OUT.(mdl) as empty table
            continue;
        end
        OUT.(mdl) = align_and_vertcat(tlist);
    end
end

% ============================== Helper: find struct ==============================
function M = extractMatchStruct(S, modelNames)
% Try common variable name 'MATCH' first; otherwise find any struct variable
% that contains at least one of the requested model fields.

    if isfield(S, 'MATCH') && isstruct(S.MATCH)
        M = S.MATCH;
        return;
    end
    % Fallback: search any struct variable with overlapping fields
    fn = fieldnames(S);
    for i = 1:numel(fn)
        val = S.(fn{i});
        if isstruct(val)
            sfields = fieldnames(val);
            if any(ismember(sfields, modelNames))
                M = val;
                return;
            end
        end
    end
    error('No suitable struct with requested model fields was found.');
end

% ============================== Helper: unique name ==============================
function out = genUniqueVarName(existing, proposal)
% Generate a variable name not in "existing" by appending numeric suffix if needed.
    out = proposal;
    idx = 2;
    while any(strcmp(existing, out))
        out = sprintf('%s_%d', proposal, idx);
        idx = idx + 1;
    end
end

% ============================== Helper: align & concat ==============================
function Tmerged = align_and_vertcat(tlist)
% Make all tables in tlist share the same variable names and classes, then vertcat.

    % Collect union of variable names and remember the first-seen class for each var
    allVars = {};
    varClass = containers.Map();
    for i = 1:numel(tlist)
        T = tlist{i};
        vnames = T.Properties.VariableNames;
        allVars = union(allVars, vnames, 'stable');
        for j = 1:numel(vnames)
            v = vnames{j};
            if ~isKey(varClass, v)
                % class of the whole column works even if T is empty
                varClass(v) = class(T.(v));
            end
        end
    end

    % Add missing variables to each table, using appropriate "missing" values
    for i = 1:numel(tlist)
        T = tlist{i};
        missingVars = setdiff(allVars, T.Properties.VariableNames, 'stable');
        if ~isempty(missingVars)
            for j = 1:numel(missingVars)
                v = missingVars{j};
                cls = varClass(v);
                T.(v) = create_missing_column(cls, height(T)); %#ok<AGROW>
            end
        end
        % Reorder columns to the union order for a clean vertcat
        T = T(:, allVars);
        tlist{i} = T; %#ok<AGROW>
    end

    % Finally, vertical concatenation
    Tmerged = vertcat(tlist{:});
end

% ============================== Helper: missing column ==============================
function col = create_missing_column(cls, m)
% Create an m-by-1 column of the specified class filled with an appropriate "missing" value.
    switch cls
        case {'double','single'}
            col = nan(m,1, cls);
        case {'logical'}
            col = false(m,1);
        case {'int8','int16','int32','int64','uint8','uint16','uint32','uint64'}
            % No NaN for integer types; use zeros
            col = zeros(m,1, cls);
        case {'string'}
            col = strings(m,1);
        case {'char'}
            % char columns inside table are rare; use blanks
            col = repmat('', m,1);
        case {'cell'}
            % Generic cell column (often used for cellstr or mixed content)
            col = cell(m,1);
            % If you expect cellstr specifically, uncomment next line:
            % col = repmat({''}, m, 1);
        case {'datetime'}
            col = NaT(m,1);
        case {'duration'}
            col = seconds(nan(m,1)); % a duration with NaN seconds
        case {'categorical'}
            col = repmat(categorical(missing), m, 1);
        otherwise
            % Fallback: create an empty cell column to be maximally permissive
            warning('create_missing_column:UnknownClass', ...
                    'Unknown variable class "%s". Using cell column as fallback.', cls);
            col = cell(m,1);
    end
end
