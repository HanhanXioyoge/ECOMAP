function out = plotRMSEBar(namesCell, rmseCell, varargin)
% plotRMSEBar
%   Plot a publication-style bar chart for RMSE values.
%
% INPUT
%   namesCell : cellstr / string array (category names)
%   rmseCell  : cell array of numeric scalars OR numeric vector (RMSE values)
%
% NAME-VALUE (optional)
%   'Title'        : figure title (default "")
%   'XLabel'       : x-axis label (default "")
%   'YLabel'       : y-axis label (default "RMSE")
%   'Sort'         : 'none'|'ascend'|'descend' (default 'none')
%   'Horizontal'   : true/false (default false)  % horizontal barh
%   'Color'        : 1x3 RGB (default Okabe–Ito blue)
%   'EdgeColor'    : 1x3 RGB (default [1 1 1])
%   'LineWidth'    : bar edge line width (default 0.8)
%   'ShowValues'   : true/false show value labels (default true)
%   'ValueFormat'  : sprintf format, e.g. '%.3f' (default '%.3g')
%   'RotateXTick'  : rotation for x tick labels (default 35)
%   'FontName'     : default 'Arial'
%   'FontSize'     : default 10
%   'YLim'         : [] or [ymin ymax] (default [])
%   'Grid'         : true/false (default false)
%
% OUTPUT (struct)
%   out.fig, out.ax, out.bar, out.names, out.rmse, out.table

% -------------------- parse options --------------------
p = inputParser;
p.addRequired('namesCell');
p.addRequired('rmseCell');

p.addParameter('Title', "", @(s)ischar(s)||isstring(s));
p.addParameter('XLabel', "", @(s)ischar(s)||isstring(s));
p.addParameter('YLabel', "RMSE", @(s)ischar(s)||isstring(s));
p.addParameter('Sort', 'none', @(s)any(strcmpi(string(s),["none","ascend","descend"])));
p.addParameter('Horizontal', false, @(b)islogical(b)&&isscalar(b));

% Okabe–Ito blue as default
p.addParameter('Color', [0.00 0.45 0.70], @(c)isnumeric(c)&&isequal(size(c),[1 3]));
p.addParameter('EdgeColor', [1 1 1], @(c)isnumeric(c)&&isequal(size(c),[1 3]));
p.addParameter('LineWidth', 0.8, @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.addParameter('ShowValues', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('ValueFormat', '%.3g', @(s)ischar(s)||isstring(s));
p.addParameter('RotateXTick', 35, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('FontName', 'Arial', @(s)ischar(s)||isstring(s));
p.addParameter('FontSize', 10, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('YLim', [], @(v)isempty(v) || (isnumeric(v)&&numel(v)==2));
p.addParameter('Grid', false, @(b)islogical(b)&&isscalar(b));

p.parse(namesCell, rmseCell, varargin{:});
opt = p.Results;

% -------------------- normalize inputs --------------------
names = string(namesCell(:));

if isnumeric(rmseCell)
    rmse = rmseCell(:);
elseif iscell(rmseCell)
    % accept {0.12, 0.3, ...} or mixed numeric
    rmse = nan(numel(rmseCell),1);
    for i = 1:numel(rmseCell)
        v = rmseCell{i};
        if isempty(v)
            rmse(i) = NaN;
        elseif isnumeric(v) && isscalar(v)
            rmse(i) = v;
        else
            rmse(i) = NaN;
        end
    end
else
    error('rmseCell must be a numeric vector or a cell array of numeric scalars.');
end

if numel(names) ~= numel(rmse)
    error('namesCell and rmse must have the same length.');
end

% drop NaN rows (optional but safer)
keep = isfinite(rmse) & strlength(names)>0;
names = names(keep);
rmse  = rmse(keep);

if isempty(names)
    error('No valid (finite) RMSE values to plot.');
end

% -------------------- sorting --------------------
switch lower(opt.Sort)
    case 'ascend'
        [rmse, ord] = sort(rmse, 'ascend');
        names = names(ord);
    case 'descend'
        [rmse, ord] = sort(rmse, 'descend');
        names = names(ord);
    case 'none'
        % keep order
end

% -------------------- plot --------------------
fig = figure('Color','w');
ax  = axes(fig);
hold(ax,'on');

% use categorical to preserve order
cats = categorical(cellstr(names));
cats = reordercats(cats, cellstr(names));

if opt.Horizontal
    b = barh(ax, cats, rmse, 0.75, 'FaceColor', opt.Color);
else
    b = bar(ax, cats, rmse, 0.75, 'FaceColor', opt.Color);
end

b.EdgeColor = opt.EdgeColor;
b.LineWidth = opt.LineWidth;

% style
ax.Color = 'w';
ax.Box = 'off';
ax.TickDir = 'out';
ax.FontName = opt.FontName;
ax.FontSize = opt.FontSize;

title(ax, string(opt.Title), 'FontName', opt.FontName, 'FontWeight','bold');
xlabel(ax, string(opt.XLabel), 'FontName', opt.FontName);
ylabel(ax, string(opt.YLabel), 'FontName', opt.FontName);

if ~opt.Horizontal
    xtickangle(ax, opt.RotateXTick);
end

if opt.Grid
    grid(ax, 'on');
else
    grid(ax, 'off');
end

% y-limits
if ~isempty(opt.YLim)
    ylim(ax, opt.YLim);
else
    if opt.Horizontal
        % xlim for barh
        xlim(ax, [0, max(rmse)*1.1]);
    else
        ylim(ax, [0, max(rmse)*1.1]);
    end
end

% value labels
if opt.ShowValues
    fmt = char(opt.ValueFormat);
    if opt.Horizontal
        % text to the right of bars
        for i = 1:numel(rmse)
            text(ax, rmse(i), cats(i), "  " + sprintf(fmt, rmse(i)), ...
                'VerticalAlignment','middle', 'HorizontalAlignment','left', ...
                'FontName', opt.FontName, 'FontSize', opt.FontSize);
        end
    else
        % text above bars
        for i = 1:numel(rmse)
            text(ax, cats(i), rmse(i), sprintf(fmt, rmse(i)), ...
                'VerticalAlignment','bottom', 'HorizontalAlignment','center', ...
                'FontName', opt.FontName, 'FontSize', opt.FontSize);
        end
    end
end

hold(ax,'off');

% -------------------- outputs --------------------
out = struct();
out.fig  = fig;
out.ax   = ax;
out.bar  = b;
out.names = names;
out.rmse  = rmse;
out.table = table(names, rmse, 'VariableNames', {'Name','RMSE'});

end
