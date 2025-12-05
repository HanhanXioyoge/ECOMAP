function [Groups, Summary, T_all_marked] = group_by_isComplex(T_all, varargin)
% GROUP_BY_ISCOMPLEX
%   Split T_all into groups by a binary column 'isComplex' (1 = Complex, 0 = NonComplex).
%   Keeps ALL columns in each group.
%
%   [Groups, Summary, T_all_marked] = group_by_isComplex(T_all)
%   [Groups, Summary, T_all_marked] = group_by_isComplex(T_all, 'Column','isComplex', 'OutDir','out')
%
% INPUT
%   T_all : table containing the binary column (default name 'isComplex').
%
% NAME-VALUE OPTIONS (optional)
%   'Column'      : char/string, the column name to classify on. Default 'isComplex'.
%   'AddLabelCol' : logical, add a categorical label column 'isComplex_label'. Default true.
%   'OutDir'      : char/string; if nonempty, write each group to files here. Default '' (no write).
%   'WriteFormat' : 'csv' | 'mat' | 'both'. Default 'csv'.
%   'FilePrefix'  : char/string, filename prefix for outputs. Default 'T_all'.
%   'Silent'      : logical, suppress console messages. Default false.
%
% OUTPUT
%   Groups       : struct with fields:
%                    - Complex     (rows where isComplex == 1/true)
%                    - NonComplex  (rows where isComplex == 0/false)
%                    - Unknown     (rows not recognized as 1/0; present only if nonempty)
%   Summary      : table with columns: Group, Count
%   T_all_marked : original table plus 'isComplex_label' (categorical) if AddLabelCol is true
%
% NOTES
%   - Robust to logical/numeric/string representations: {1,true,'1','true','yes'} => Complex;
%     {0,false,'0','false','no'} => NonComplex; others => Unknown.
%
% Author: (your name), 2025-11-07

    % ------------------ Parse inputs ------------------
    ip = inputParser;
    ip.addParameter('Column',      'isComplex', @(x)ischar(x)||isstring(x));
    ip.addParameter('AddLabelCol', true,        @(x)islogical(x)&&isscalar(x));
    ip.addParameter('OutDir',      '',          @(x)ischar(x)||isstring(x));
    ip.addParameter('WriteFormat', 'csv',       @(x)ischar(x)||isstring(x));
    ip.addParameter('FilePrefix',  'T_all',     @(x)ischar(x)||isstring(x));
    ip.addParameter('Silent',      false,       @(x)islogical(x)&&isscalar(x));
    ip.parse(varargin{:});

    colName    = char(ip.Results.Column);
    addLabel   = ip.Results.AddLabelCol;
    outDir     = char(ip.Results.OutDir);
    fmt        = lower(char(ip.Results.WriteFormat));
    filePrefix = char(ip.Results.FilePrefix);
    silent     = ip.Results.Silent;

    if ~istable(T_all)
        error('T_all must be a table.');
    end
    if ~ismember(colName, T_all.Properties.VariableNames)
        error('Column "%s" not found in T_all.', colName);
    end

    % ------------------ Normalize to labels ------------------
    raw = T_all.(colName);
    labels = classify_iscomplex(raw);  % string array: "Complex" | "NonComplex" | "Unknown"

    % ------------------ Add label column (optional) ------------------
    if addLabel
        T_all_marked = add_or_replace_var(T_all, 'isComplex_label', categorical(labels));
    else
        T_all_marked = T_all;
    end

    % ------------------ Split by labels ------------------
    maskC = labels == "Complex";
    maskN = labels == "NonComplex";
    maskU = labels == "Unknown";

    Groups = struct();
    Groups.Complex    = T_all(maskC, :);
    Groups.NonComplex = T_all(maskN, :);
    if any(maskU)
        Groups.Unknown = T_all(maskU, :);
    end

    % ------------------ Summary table ------------------
    grpNames = ["Complex","NonComplex"];
    counts   = [sum(maskC), sum(maskN)];
    if any(maskU)
        grpNames = [grpNames, "Unknown"];
        counts   = [counts,   sum(maskU)];
    end
    Summary = table(categorical(grpNames.', grpNames), counts.', ...
                    'VariableNames', {'Group','Count'});

    % ------------------ Optional write ------------------
    if ~isempty(outDir)
        if ~exist(outDir, 'dir'), mkdir(outDir); end
        write_one(Groups.Complex,    fullfile(outDir, [filePrefix '_complex']),    fmt, silent);
        write_one(Groups.NonComplex, fullfile(outDir, [filePrefix '_noncomplex']), fmt, silent);
        if isfield(Groups, 'Unknown')
            write_one(Groups.Unknown, fullfile(outDir, [filePrefix '_unknown']), fmt, silent);
        end
    end
end

% ===================== helpers =====================

function labels = classify_iscomplex(col)
% Map various representations to "Complex" / "NonComplex" / "Unknown".
    n = numel(col);
    labels = strings(n,1);

    for i = 1:n
        v = col(i);
        lab = "Unknown";
        if islogical(v)
            lab = tern(v, "Complex", "NonComplex");
        elseif isnumeric(v)
            if isfinite(v) && ~isnan(v)
                if v == 1
                    lab = "Complex";
                elseif v == 0
                    lab = "NonComplex";
                end
            end
        elseif isstring(v) || ischar(v)
            s = lower(strtrim(string(v)));
            if any(strcmp(s, ["1","true","yes","y","complex"]))
                lab = "Complex";
            elseif any(strcmp(s, ["0","false","no","n","noncomplex","non-complex"]))
                lab = "NonComplex";
            end
        elseif iscell(v) && ~isempty(v)
            % Handle cell content recursively (rare)
            lab = classify_iscomplex(v{1});
        end
        labels(i) = lab;
    end
end

function out = tern(cond, a, b)
% Ternary helper
    if cond, out = a; else, out = b; end
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

function write_one(T, basePath, fmt, silent)
% Write table T as csv/mat/both using basePath (no extension)
    switch fmt
        case 'csv'
            writetable(T, [basePath '.csv']);
            if ~silent, fprintf('[group_by_isComplex] Wrote %s.csv (%d rows)\n', basePath, height(T)); end
        case 'mat'
            S = T; %#ok<NASGU>
            save([basePath '.mat'], 'S');
            if ~silent, fprintf('[group_by_isComplex] Wrote %s.mat (%d rows)\n', basePath, height(T)); end
        case 'both'
            writetable(T, [basePath '.csv']);
            S = T; %#ok<NASGU>
            save([basePath '.mat'], 'S');
            if ~silent, fprintf('[group_by_isComplex] Wrote %s.[csv|mat] (%d rows)\n', basePath, height(T)); end
        otherwise
            error('Unknown WriteFormat "%s". Use csv|mat|both.', fmt);
    end
end
