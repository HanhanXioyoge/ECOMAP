function T = rxnStructToTable(S, varargin)
%RXNSTRUCTTOTABLE Convert a struct into a MATLAB table for display.
%   T = rxnStructToTable(S) expects S to contain (either as fields or per
%   element if S is a struct array) the items:
%       - rxns        : reaction ids/names
%       - substrates  : each row can be a string/char, or a cell of strings
%       - eccodes     : EC numbers; single or multiple per row are ok
%       - kcats       : numeric vector OR per-row values (numeric/string/cell)
%
%   Name-Value options:
%       'AsGUI'      : true/false (default false). If true, opens a figure
%                      with a uitable to view the result.
%       'FigureName' : title for the GUI window (default 'Reaction Table')
%       'Delimiter'  : joiner for multi-valued cells (default '; ')
%
%   Returns a MATLAB table T with variables: rxns, substrates, eccodes, kcats.

% ------- parse options -------
p = inputParser;
p.addParameter('AsGUI', false, @(x)islogical(x) && isscalar(x));
p.addParameter('FigureName','Reaction Table', @(x)ischar(x)||isstring(x));
p.addParameter('Delimiter','; ', @(x)ischar(x)||isstring(x));
p.parse(varargin{:});
opt = p.Results;
joiner = char(opt.Delimiter);

% ------- pick columns from struct or struct array -------
rxns_raw       = pickField(S, 'rxns');
subs_raw       = pickField(S, 'substrates');
ec_raw         = pickField(S, 'eccodes');
kcats_raw      = pickField(S, 'kcats');

% Determine row count
n = max([len(rxns_raw), len(subs_raw), len(ec_raw), len(kcats_raw)]);
if n == 0
    T = table();
    warning('No data found in fields rxns/substrates/eccodes/kcats.');
    return;
end

% ------- normalize each column to column vectors -------
rxns_col  = normalizeColumn(rxns_raw, n, joiner);
subs_col  = normalizeColumn(subs_raw, n, joiner);
ec_col    = normalizeColumn(ec_raw,   n, joiner);

% For kcats, try to keep numeric if it is a numeric vector of length n
if isnumeric(kcats_raw) && isvector(kcats_raw) && numel(kcats_raw)==n
    kcats_col = reshape(kcats_raw, [], 1);
else
    kcats_col = normalizeColumn(kcats_raw, n, joiner);
end

% ------- build table -------
T = table(rxns_col, subs_col, ec_col, kcats_col, ...
    'VariableNames', {'ReactionName','Substrate','eccodes','kcat_value'});

% nicer string type for display (keep kcats原样；数值会保持为double，字符串保持字符串)
T.ReactionName      = toStringColumn(T.ReactionName);
T.Substrate         = toStringColumn(T.Substrate);
T.eccodes           = toStringColumn(T.eccodes);

% ------- NEW: filter out rows with kcat == 0 -------
T = filterZeroKcat(T);

% ------- optional GUI -------
if opt.AsGUI
    f = figure('Name', opt.FigureName, 'NumberTitle','off', ...
               'Color','w', 'Position',[100 100 900 520]);
    uitable(f, 'Data', table2cell(T), ...
               'ColumnName', T.Properties.VariableNames, ...
               'Units','normalized','Position',[0 0 1 1]);
end
end

% ================= helpers =================
function v = pickField(S, field)
if ~isstruct(S)
    v = [];
    return;
end
if isscalar(S)
    if isfield(S, field), v = S.(field); else, v = []; end
else
    hasF = arrayfun(@(x)isfield(x,field), S);
    if ~all(hasF)
        v = cell(numel(S),1); v(:) = {''};
        return;
    end
    try
        v = arrayfun(@(x)x.(field), S, 'UniformOutput', false);
    catch
        v = arrayfun(@(x)tryget(x,field), S, 'UniformOutput', false);
    end
end
end

function out = tryget(s, f)
try
    out = s.(f);
catch
    out = '';
end
end

function n = len(x)
if isempty(x), n = 0; return; end
if iscell(x) || isstring(x)
    n = numel(x);
elseif ischar(x)
    if isrow(x), n = 1; else, n = size(x,1); end
elseif isnumeric(x) || islogical(x)
    n = numel(x);
else
    n = numel(x);
end
end

function c = normalizeColumn(x, n, joiner)
toStr = @(s) elemToStr(s, joiner);
if isempty(x)
    c = repmat({''}, n, 1);
    return;
end
if isstring(x)
    x = cellstr(x(:));
elseif ischar(x)
    if isrow(x), x = {x}; else, x = cellstr(x); end
elseif isnumeric(x) || islogical(x)
    x = arrayfun(@(v) num2str(v), x(:), 'UniformOutput', false);
elseif iscell(x)
    x = x(:);
else
    x = cellstr(string(x(:)));
end
c = cellfun(toStr, x, 'UniformOutput', false);
m = numel(c);
if m < n, c = [c; repmat({''}, n-m, 1)]; end
if m > n, c = c(1:n); end
end

function s = elemToStr(v, joiner)
if isstring(v)
    s = char(v);
elseif ischar(v)
    s = v;
elseif isnumeric(v) || islogical(v)
    if isscalar(v), s = num2str(v);
    else, s = ['[', strjoin(arrayfun(@num2str, v(:)','UniformOutput',false), ','), ']'];
    end
elseif iscell(v)
    s = strjoin(cellfun(@(w) elemToStr(w, joiner), v, 'UniformOutput', false), joiner);
else
    s = char(string(v));
end
s = s(:).';
end

function col = toStringColumn(col)
if iscell(col)
    col = string(col);
elseif ischar(col)
    col = string(col);
end
end

function T = filterZeroKcat(T)
% Remove rows whose kcats is numerically 0.
% Works for numeric or string kcats.
if isempty(T), return; end
if isnumeric(T.kcat_value)
    keep = T.kcat_value ~= 0;
else
    knum = str2double(string(T.kcat_value)); % non-numeric => NaN
    keep = ~(knum == 0);
    keep(isnan(knum)) = true;          % keep non-parsable (treat as not-zero)
end
T = T(keep, :);
end
