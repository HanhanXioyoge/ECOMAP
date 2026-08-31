function result = mdpBuildOkoIntervals(predictions, predictor)
%MDPBUILDOKOINTERVALS  Build OKO+ kcat intervals from DL predictor tables.
%
%   result = mdpBuildOkoIntervals(predictions)
%   result = mdpBuildOkoIntervals(predictions, predictor)
%
%   Inputs:
%     predictions -- cell array of structs. Each struct must expose ``table``
%                    (a MATLAB table) OR a pair of fields ``columns`` (cellstr)
%                    and ``rows`` (cell array of cellstr/double rows).
%     predictor   -- optional char filter; only rows with matching predictor
%                    name are kept. Empty (default) keeps all rows.
%
%   Returns the standard bridge envelope (see CONTRACT.md) carrying:
%     .columns / .rows -- the aggregated interval table.
%
%   Defers to buildOkoIntervals in scripts/external_kcat_prediction/.
    here = fileparts(mfilename('fullpath'));
    addpath_once(here);
    addpath_once(fullfile(here, '..', '..', 'external_kcat_prediction'));

    if nargin < 1 || isempty(predictions)
        result = make_err('err_param_invalid', 'predictions cell array required');
        return;
    end
    if nargin < 2, predictor = ''; end

    tables = cell(1, numel(predictions));
    for i = 1:numel(predictions)
        entry = predictions{i};
        if istable(entry)
            tables{i} = entry;
        elseif isstruct(entry) && isfield(entry, 'rows') && isfield(entry, 'columns')
            tables{i} = cellTableFromColumns(entry.columns, entry.rows);
        else
            result = make_err('err_param_invalid', ...
                sprintf('prediction %d must be a table or {rows,columns}', i));
            return;
        end
    end

    args = {tables};
    if ~isempty(predictor)
        args = [args, {'Predictor', char(predictor)}];
    end

    try
        if exist('buildOkoIntervals', 'file') ~= 2 && exist('buildOkoIntervals', 'file') ~= 5
            error('buildOkoIntervals:notFound', 'buildOkoIntervals helper not on path');
        end
        intervals = buildOkoIntervals(args{:});
    catch err
        result = make_err('err_kcat_merge', err.message);
        return;
    end

    [cols, rows] = tableToColumns(intervals);
    result = make_ok(struct('columns', {cols}, 'rows', {rows}));
end

function T = cellTableFromColumns(cols, rows)
    nRows = numel(rows);
    varNames = cell(1, numel(cols));
    varCells = cell(1, numel(cols));
    for j = 1:numel(cols)
        col = cols{j};
        if ischar(col) || isstring(col)
            varNames{j} = char(col);
        else
            varNames{j} = sprintf('col%d', j);
        end
        vals = cell(nRows, 1);
        for i = 1:nRows
            v = rows{i};
            if numel(v) >= j
                vals{i} = v{j};
            else
                vals{i} = [];
            end
        end
        varCells{j} = vals;
    end
    T = table(varCells{:}, 'VariableNames', varNames);
end

function [cols, rows] = tableToColumns(T)
    if isempty(T)
        cols = {};
        rows = {};
        return;
    end
    cols = T.Properties.VariableNames;
    cells = table2cell(T);
    n = height(T);
    out = cell(1, n);
    for i = 1:n, out{i} = cells(i, :); end
    rows = out;
end