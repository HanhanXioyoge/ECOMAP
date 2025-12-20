function out = plotMultiECDF_Log10(dataGroups, labels, varargin)
% plotMultiECDF_Log10
%   Plot ECDFs for multiple datasets in ONE figure with log10 x-axis
%   (ticks shown as 10^n), a y=0.5 horizontal line, and for each group:
%   a vertical dashed line from y=0 to y=0.5 at the intersection point.
%   Also draw a full box frame WITHOUT showing top/right ticks by overlaying
%   a second transparent axes.
%
% INPUT
%   dataGroups : cell array {x1, x2, ...} OR numeric matrix (each column one group)
%   labels     : cellstr/string array, same number as groups
%
% NAME-VALUE (optional)
%   'Colors'   : colors for each group
%               - Nx3 RGB matrix (each row one color)
%               - 1x3 RGB row (reuse for all)
%               - cell array like {'r',[0 0.5 0.2],'#1f77b4'} (cycled if shorter)
%               - [] or omitted: use built-in Okabe–Ito palette (fixed order)
%   'Title'    : figure title (default '')
%   'XLabel'   : x-axis label (default 'Value')
%   'YLabel'   : y-axis label (default 'Empirical CDF')
%   'ShowGrid' : true/false (default false, i.e., grid hidden)
%   'LineWidth': line width (default 2)
%   'ShowHalf' : true/false show y=0.5 line (default true)
%
% OUTPUT (struct)
%   out.fig, out.ax, out.axBox, out.lines, out.vlines, out.usedData, out.xAtHalf

% -------------------- parse inputs --------------------
p = inputParser;
p.addRequired('dataGroups');
p.addRequired('labels');
p.addParameter('Colors', [], @(c) isnumeric(c) || iscell(c) || isstring(c) || ischar(c));
p.addParameter('Title', '', @(s) isstring(s) || ischar(s));
p.addParameter('XLabel', 'Value', @(s) isstring(s) || ischar(s));
p.addParameter('YLabel', 'Empirical CDF', @(s) isstring(s) || ischar(s));
p.addParameter('ShowGrid', false, @(x) islogical(x) && isscalar(x));
p.addParameter('LineWidth', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('ShowHalf', true, @(x) islogical(x) && isscalar(x));
p.parse(dataGroups, labels, varargin{:});
opt = p.Results;

% -------------------- normalize groups --------------------
if isnumeric(dataGroups)
    G = size(dataGroups, 2);
    tmp = cell(1, G);
    for i = 1:G
        tmp{i} = dataGroups(:, i);
    end
    dataGroups = tmp;
elseif ~iscell(dataGroups)
    error('dataGroups must be a cell array or a numeric matrix (columns = groups).');
end

labels = string(labels);
nG = numel(dataGroups);
if numel(labels) ~= nG
    error('Number of labels must match number of data groups.');
end

% -------------------- create figure --------------------
fig = figure('Color', 'w');
ax  = axes(fig); hold(ax, 'on');

% Built-in palette (Okabe–Ito), fixed order, colorblind-friendly
okabeIto = [
    0.9020 0.6235 0.0000  % orange
    0.3373 0.7059 0.9137  % sky blue
    0.0000 0.6196 0.4510  % bluish green
    0.9412 0.8941 0.2588  % yellow
    0.0000 0.4471 0.6980  % blue
    0.8353 0.3686 0.0000  % vermillion
    0.8000 0.4745 0.6549  % reddish purple
    0.0000 0.0000 0.0000  % black
];

% Prepare colors (cycle if shorter)
colorList = localNormalizeColors(opt.Colors, nG, okabeIto);

usedData = cell(nG, 1);
xAtHalf  = nan(nG, 1);

lines  = gobjects(nG, 1);
vlines = gobjects(nG, 1);
marks  = gobjects(nG, 1);
xmin = inf; xmax = -inf;

% -------------------- plot each ECDF --------------------
for i = 1:nG
    x = dataGroups{i};
    x = x(:);
    x = x(isfinite(x));      % drop NaN/Inf
    x = x(x > 0);            % log-x requires positive values

    if isempty(x)
        warning('Group %d (%s) has no positive finite values. Skipped.', i, labels(i));
        continue;
    end

    usedData{i} = x;

    [f, xx] = ecdf(x);

    % ECDF curve
    lines(i) = stairs(ax, xx, f, ...
        'LineWidth', opt.LineWidth, ...
        'DisplayName', labels(i), ...
        'Color', colorList{i});

    xmin = min(xmin, min(xx));
    xmax = max(xmax, max(xx));

    % x at y=0.5 (intersection)
    xHalf = interp1(f, xx, 0.5, 'linear', 'extrap');
    xAtHalf(i) = xHalf;

    marks(i) = plot(ax, xHalf, 0.5, 'o', ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', 'none', ... 
        'MarkerEdgeColor', colorList{i}, ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

    % vertical dashed line down to x-axis (y=0) from y=0.5
    vlines(i) = line(ax, [xHalf xHalf], [0 0.5], ...
        'LineStyle', '--', ...
        'LineWidth', max(1, opt.LineWidth-0.5), ...
        'Color', colorList{i}, ...
        'HandleVisibility', 'off');
end

% -------------------- axis formatting --------------------
set(ax, 'XScale', 'log');

if isfinite(xmin) && isfinite(xmax) && xmin > 0 && xmax > 0
    emin = floor(log10(xmin));
    emax = ceil(log10(xmax));
    exps = emin:emax;

    if numel(exps) > 8
        step = ceil(numel(exps)/8);
        exps = exps(1:step:end);
    end

    ax.XTick = 10.^exps;       % ticks at powers of 10 (shown as 10^n)
    xlim(ax, [10^emin, 10^emax]);
end

ylim(ax, [0, 1]);
xlabel(ax, opt.XLabel);
ylabel(ax, opt.YLabel);
title(ax, opt.Title);

% grid default hidden
if opt.ShowGrid
    grid(ax, 'on');
else
    grid(ax, 'off');
end

% y=0.5 line
if opt.ShowHalf
    yline(ax, 0.5, ':', '0.5', 'HandleVisibility', 'off');
end

legend(ax, 'Location', 'best', 'Interpreter', 'none');

% -------------------- box without top/right ticks (overlay axes) --------------------
ax.Box = 'off';  % keep only bottom/left ticks on main axes
ax.TickDir = 'out';

axBox = axes('Parent', fig, ...
    'Position', ax.Position, ...
    'Color', 'none', ...
    'Box', 'on', ...
    'XTick', [], 'YTick', [], ...          % remove all ticks on overlay axes
    'XLim', ax.XLim, 'YLim', ax.YLim, ...
    'XScale', ax.XScale, 'YScale', ax.YScale, ...
    'HitTest', 'off', 'HandleVisibility', 'off');

linkaxes([ax axBox], 'xy');
uistack(ax, 'top');  % ensure main axes stays on top for labels/legend

% -------------------- outputs --------------------
out = struct();
out.fig      = fig;
out.ax       = ax;
out.axBox    = axBox;
out.lines    = lines;
out.vlines   = vlines;
out.usedData = usedData;
out.xAtHalf  = xAtHalf;

end

% ===== helper =====
function colorList = localNormalizeColors(colorsIn, nG, defaultPalette)
% returns cell array {c1,c2,...} each c is RGB 1x3
    if isempty(colorsIn)
        % Use built-in Okabe–Ito (fixed order), cycle if needed
        colorList = cell(nG, 1);
        for i = 1:nG
            colorList{i} = defaultPalette(mod(i-1, size(defaultPalette,1)) + 1, :);
        end
        return
    end

    % If user gave a char/string (single color name or hex), wrap it
    if ischar(colorsIn) || isstring(colorsIn)
        colorsIn = {colorsIn};
    end

    if isnumeric(colorsIn)
        if size(colorsIn,2) ~= 3
            error('Colors as numeric must be Nx3 RGB.');
        end
        if size(colorsIn,1) == 1
            colorsIn = repmat(colorsIn, nG, 1);
        end
        colorList = cell(nG,1);
        for i = 1:nG
            colorList{i} = colorsIn(mod(i-1,size(colorsIn,1))+1, :);
        end
        return
    end

    if iscell(colorsIn)
        colorList = cell(nG,1);
        for i = 1:nG
            ci = colorsIn{mod(i-1,numel(colorsIn))+1};
            colorList{i} = localToRGB(ci);
        end
        return
    end

    error('Unsupported Colors format.');
end

function rgb = localToRGB(ci)
    if isnumeric(ci)
        if ~isequal(size(ci), [1 3])
            error('RGB color must be 1x3.');
        end
        rgb = ci;
        return
    end

    ci = string(ci);
    if startsWith(ci, "#")
        rgb = sscanf(extractAfter(ci,1), '%2x%2x%2x', [1 3]) / 255;
        return
    end

    % short names: r g b c m y k w
    switch lower(char(ci))
        case 'r', rgb = [1 0 0];
        case 'g', rgb = [0 1 0];
        case 'b', rgb = [0 0 1];
        case 'c', rgb = [0 1 1];
        case 'm', rgb = [1 0 1];
        case 'y', rgb = [1 1 0];
        case 'k', rgb = [0 0 0];
        case 'w', rgb = [1 1 1];
        otherwise
            error('Unrecognized color spec: %s', ci);
    end
end
