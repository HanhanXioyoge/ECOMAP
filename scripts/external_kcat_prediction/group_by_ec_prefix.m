function [Groups, Summary, T_all_marked] = group_by_ec_prefix(T_all, varargin)
% GROUP_BY_EC_PREFIX
%   Group rows of T_all by the first EC class digit in column 'ec'.
%   Example: "2.7.1.35" -> class "2", grouped under field "EC2".
%
%   [Groups, Summary, T_all_marked] = group_by_ec_prefix(T_all)
%   [Groups, Summary, T_all_marked] = group_by_ec_prefix(T_all, 'ECColumn','ec', 'OutDir','results')
%
% INPUT
%   T_all : table, must contain an EC column (default name 'ec').
%
% NAME-VALUE OPTIONS (all optional)
%   'ECColumn'     : char/string, EC column name. Default 'ec'.
%   'AddClassCol'  : logical, add 'EC_Class' column to the returned table. Default true.
%   'OutDir'       : char/string; if nonempty, write each group to files in this folder. Default '' (no write).
%   'WriteFormat'  : 'csv' | 'mat' | 'both'. Default 'csv'.
%   'FilePrefix'   : char/string, prefix for output files. Default 'T_all'.
%   'Silent'       : logical, suppress console messages. Default false.
%
% OUTPUT
%   Groups       : struct with fields like EC1..EC7 and EC_Unknown (only present if nonempty),
%                  each field is a table containing ALL columns from T_all for that class.
%   Summary      : table with columns: EC_Class, Count; sorted in EC order 1..7 then Unknown.
%   T_all_marked : T_all with an added 'EC_Class' column (categorical) if 'AddClassCol'==true.
%
% NOTES
%   - Robust to EC values with prefixes ('EC 2.7.1.35'), mixed case, leading/trailing spaces,
%     or multiple ECs separated by punctuation; the function uses the first leading digit it finds.
%   - Rows without a detectable leading digit are classified as "EC_Unknown".
%
% Author: (your name), 2025-11-07

    % ------------------ Parse inputs ------------------
    ip = inputParser;
    ip.addParameter('ECColumn',    'ec', @(x)ischar(x)||isstring(x));
    ip.addParameter('AddClassCol', true, @(x)islogical(x)&&isscalar(x));
    ip.addParameter('OutDir',      '',   @(x)ischar(x)||isstring(x));
    ip.addParameter('WriteFormat', 'csv', @(x)ischar(x)||isstring(x));
    ip.addParameter('FilePrefix',  'T_all', @(x)ischar(x)||isstring(x));
    ip.addParameter('Silent',      false, @(x)islogical(x)&&isscalar(x));
    ip.parse(varargin{:});

    ecColName  = char(ip.Results.ECColumn);
    addClass   = ip.Results.AddClassCol;
    outDir     = char(ip.Results.OutDir);
    fmt        = lower(char(ip.Results.WriteFormat));
    filePrefix = char(ip.Results.FilePrefix);
    silent     = ip.Results.Silent;

    if ~istable(T_all)
        error('T_all must be a table.');
    end
    if ~ismember(ecColName, T_all.Properties.VariableNames)
        error('EC column "%s" not found in T_all.', ecColName);
    end

    % ------------------ Normalize EC strings & extract class ------------------
    ECstr = to_str_col(T_all.(ecColName));
    ECstr = strip(regexprep(ECstr, '^\s*EC\s*', '', 'ignorecase')); % drop leading 'EC'
    % If multiple ECs present, we still only need the first leading digit
    % Extract the very first digit run at the string start (after trimming)
    firstDigit = strings(size(ECstr));
    for i = 1:numel(ECstr)
        s = ECstr(i);
        % Trim to first token if separated by comma/semicolon/pipe/space
        % (this helps when someone writes "2.7.1.35; 1.1.1.1")
        s = string(regexp(s, '^[^,;| ]+', 'match', 'once'));
        m = regexp(s, '^\s*(\d+)', 'tokens', 'once'); % first leading digits
        if ~isempty(m)
            firstDigit(i) = string(m{1});
        else
            % As a fallback, look for a digit anywhere early in the string
            m2 = regexp(s, '(\d+)', 'tokens', 'once');
            if ~isempty(m2)
                firstDigit(i) = string(m2{1});
            else
                firstDigit(i) = "";
            end
        end
    end

    % Map to EC class labels
    EC_Class = strings(size(firstDigit));
    for i = 1:numel(firstDigit)
        d = str2double(firstDigit(i));
        if isfinite(d) && ismember(d, 1:7)
            EC_Class(i) = "EC" + string(d);
        else
            EC_Class(i) = "EC_Unknown";
        end
    end

    % ------------------ Attach class column (optional) ------------------
    if addClass
        T_all_marked = add_or_replace_var(T_all, 'EC_Class', categorical(EC_Class));
    else
        T_all_marked = T_all;
    end

    % ------------------ Split into groups, keep ALL columns ------------------
    catsDesired = ["EC1","EC2","EC3","EC4","EC5","EC6","EC7","EC_Unknown"];
    presentCats = unique(EC_Class, 'stable');
    presentCats = presentCats(~ismissing(presentCats)); % drop <missing>

    Groups = struct();
    for c = 1:numel(presentCats)
        cname = presentCats(c);
        mask  = EC_Class == cname;
        field = matlab.lang.makeValidName(char(cname));
        Groups.(field) = T_all(mask, :); %#ok<STRNU>
    end

    % ------------------ Build summary table ------------------
    allCounts = table();
    for c = 1:numel(catsDesired)
        cname = catsDesired(c);
        cnt = sum(EC_Class == cname);
        if cnt > 0
            allCounts = [allCounts; table(categorical(cname), cnt, 'VariableNames', {'EC_Class','Count'})]; %#ok<AGROW>
        end
    end
    Summary = allCounts;

    % ------------------ Optional writing to disk ------------------
    if ~isempty(outDir)
        if ~exist(outDir, 'dir'), mkdir(outDir); end
        fns = fieldnames(Groups);
        for i = 1:numel(fns)
            fname = sprintf('%s_%s', filePrefix, fns{i});
            switch fmt
                case 'csv'
                    writetable(Groups.(fns{i}), fullfile(outDir, [fname '.csv']));
                case 'mat'
                    S = Groups.(fns{i}); %#ok<NASGU>
                    save(fullfile(outDir, [fname '.mat']), 'S');
                case 'both'
                    writetable(Groups.(fns{i}), fullfile(outDir, [fname '.csv']));
                    S = Groups.(fns{i}); %#ok<NASGU>
                    save(fullfile(outDir, [fname '.mat']), 'S');
                otherwise
                    error('Unknown WriteFormat "%s". Use csv|mat|both.', fmt);
            end
            if ~silent
                fprintf('[group_by_ec_prefix] Wrote %s (%d rows)\n', fns{i}, height(Groups.(fns{i})));
            end
        end
    end
end

% ===================== helpers =====================

function s = to_str_col(col)
% Convert a table column to string array robustly.
    if isstring(col)
        s = col;
    elseif iscellstr(col)
        s = string(col);
    elseif iscell(col)
        s = strings(size(col));
        for i = 1:numel(col)
            if isstring(col{i})
                s(i) = col{i};
            elseif ischar(col{i})
                s(i) = string(col{i});
            elseif isnumeric(col{i}) || islogical(col{i})
                s(i) = string(col{i});
            else
                s(i) = "";
            end
        end
    elseif isnumeric(col) || islogical(col)
        s = string(col);
    else
        s = string(col);
    end
end

function T2 = add_or_replace_var(T, varName, col)
% Add new variable or replace if it already exists.
    if ismember(varName, T.Properties.VariableNames)
        T.(varName) = col;
        T2 = T;
    else
        T2 = addvars(T, col, 'NewVariableNames', varName, 'After', 1); % insert near the front
    end
end
