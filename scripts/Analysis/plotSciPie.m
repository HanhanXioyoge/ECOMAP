function [fig, ax, stats] = plotSciPie(x, varargin)
% plotSciPie  Scientific-style pie chart with FIXED colors for known types.
%
% Example:
%   x = {'brenda','median','CatPred','custom','CompleteMatch','brenda'};
%   plotSciPie(x, 'Title','kcat source', 'Order','descend');
%
% Inputs:
%   x : cell / string / categorical / numeric (converted to string)
%
% Name-Value Options:
%   'Title'             : char/string, default ""
%   'CaseSensitive'     : true/false, default true
%   'Order'             : 'descend'|'ascend'|'stable', default 'descend'
%   'OtherThreshold'    : merge categories with Percent < threshold into "Other" (0 disables), default 0
%   'TopN'              : keep only top N categories; rest merged into "Other" (0 disables), default 0
%   'LabelMode'         : 'percent'|'name+percent'|'name+percent+count', default 'name+percent'
%   'MinPercentForLabel': hide labels smaller than this percent, default 1
%   'ShowLegend'        : true/false, default true
%   'ExplodeLargest'    : true/false, default false
%   'HideAxes'          : true/false, default true
%   'FontName'          : default 'Arial'
%   'FontSize'          : default 10
%
% Outputs:
%   fig   : figure handle
%   ax    : axes handle
%   stats : table with Category, Count, Percent

% -------------------- Parse options --------------------
p = inputParser;
p.addParameter('Title', "", @(s)ischar(s)||isstring(s));
p.addParameter('CaseSensitive', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('Order', 'descend', @(s)any(strcmpi(string(s),["descend","ascend","stable"])));
p.addParameter('OtherThreshold', 0, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('TopN', 0, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&mod(v,1)==0);
p.addParameter('LabelMode', 'name+percent', @(s)any(strcmpi(string(s),["percent","name+percent","name+percent+count"])));
p.addParameter('MinPercentForLabel', 1, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('ShowLegend', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('ExplodeLargest', false, @(b)islogical(b)&&isscalar(b));
p.addParameter('HideAxes', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('FontName', 'Arial', @(s)ischar(s)||isstring(s));
p.addParameter('FontSize', 10, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.parse(varargin{:});
opt = p.Results;

% -------------------- Sanitize input to string --------------------
if iscategorical(x)
    s_raw = string(x(:));
elseif isstring(x)
    s_raw = x(:);
elseif iscell(x)
    try
        s_raw = string(x(:));
    catch
        % fallback for mixed cell types
        s_raw = strings(numel(x),1);
        for i = 1:numel(x)
            if isempty(x{i})
                s_raw(i) = "";
            elseif ischar(x{i}) || isstring(x{i}) || isnumeric(x{i}) || islogical(x{i})
                s_raw(i) = string(x{i});
            else
                s_raw(i) = string(class(x{i}));
            end
        end
    end
else
    s_raw = string(x(:));
end

% remove missing/empty
s_raw = strip(s_raw);
s_raw = s_raw(~ismissing(s_raw) & strlength(s_raw)>0);
if isempty(s_raw)
    error('plotSciPie:EmptyData', 'Input contains no valid category values.');
end

% -------------------- Count categories (preserve display label) --------------------
if opt.CaseSensitive
    key_raw = s_raw;               % keys used for counting
    [keys, ia, idx] = unique(key_raw, 'stable');
    labels = keys;                 % display labels
else
    key_raw = lower(s_raw);
    [keys, ia, idx] = unique(key_raw, 'stable'); % keys are lowercased
    labels = s_raw(ia);                              % display label = first observed original
end
counts = accumarray(idx, 1, [numel(keys), 1]);

% order
switch lower(opt.Order)
    case 'descend'
        [counts, ord] = sort(counts, 'descend');
        keys   = keys(ord);
        labels = labels(ord);
    case 'ascend'
        [counts, ord] = sort(counts, 'ascend');
        keys   = keys(ord);
        labels = labels(ord);
    case 'stable'
        % keep as-is
end

% -------------------- Merge to "Other" if needed --------------------
% helper to append/merge "Other"
    function [keys2, labels2, counts2] = mergeOther(keys1, labels1, counts1, maskToOther)
        if ~any(maskToOther), keys2=keys1; labels2=labels1; counts2=counts1; return; end
        otherCount = sum(counts1(maskToOther));
        keepMask = ~maskToOther;
        keys2   = [keys1(keepMask); "other"];
        labels2 = [labels1(keepMask); "Other"];
        counts2 = [counts1(keepMask); otherCount];
    end

% TopN merge
if opt.TopN > 0 && opt.TopN < numel(counts)
    maskToOther = true(size(counts));
    maskToOther(1:opt.TopN) = false;
    [keys, labels, counts] = mergeOther(keys, labels, counts, maskToOther);
end

% Percent threshold merge
total = sum(counts);
perc  = counts / total * 100;
if opt.OtherThreshold > 0
    maskToOther = perc < opt.OtherThreshold;
    % keep at least one category besides Other
    if any(maskToOther) && sum(~maskToOther) >= 1
        [keys, labels, counts] = mergeOther(keys, labels, counts, maskToOther);
        total = sum(counts);
        perc  = counts / total * 100;
    end
end

% -------------------- Build slice labels --------------------
sliceText = strings(size(labels));
for i = 1:numel(labels)
    if perc(i) < opt.MinPercentForLabel
        sliceText(i) = "";
        continue;
    end
    switch lower(opt.LabelMode)
        case 'percent'
            sliceText(i) = sprintf('%.1f%%', perc(i));
        case 'name+percent'
            sliceText(i) = sprintf('%s (%.1f%%)', labels(i), perc(i));
        case 'name+percent+count'
            sliceText(i) = sprintf('%s (%.1f%%, n=%d)', labels(i), perc(i), counts(i));
    end
end

% -------------------- FIXED color map for known types --------------------
% Use lower-case keys for lookup (so "CatPred" and "catpred" share color)
colorMap = containers.Map('KeyType','char','ValueType','any');

% Okabe–Ito palette (colorblind-friendly), assigned to your fixed types:
colorMap('brenda')        = [0.00 0.45 0.70]; % blue
colorMap('median')        = [0.90 0.62 0.00]; % orange
colorMap('catpred')       = [0.00 0.62 0.45]; % green/teal
colorMap('custom')        = [0.84 0.37 0.00]; % vermillion
colorMap('completematch') = [0.80 0.47 0.65]; % purple
colorMap('other')         = [0.60 0.60 0.60]; % grey for merged small classes

% fallback colors for unknown categories (if any appear)
fallback = [
    0.35 0.70 0.90
    0.94 0.89 0.26
    0.20 0.20 0.20
    0.57 0.57 0.57
    0.12 0.47 0.71
    1.00 0.50 0.05
    0.17 0.63 0.17
    0.58 0.40 0.74
];
fallbackIdx = 1;

sliceColors = zeros(numel(keys), 3);
for i = 1:numel(keys)
    k = lower(char(keys(i)));
    if isKey(colorMap, k)
        sliceColors(i,:) = colorMap(k);
    else
        sliceColors(i,:) = fallback(fallbackIdx, :);
        fallbackIdx = fallbackIdx + 1;
        if fallbackIdx > size(fallback,1), fallbackIdx = 1; end
    end
end

% -------------------- Plot --------------------
fig = figure('Color','w');
ax = axes(fig);
hold(ax, 'on');

explode = zeros(size(counts));
if opt.ExplodeLargest && ~isempty(explode)
    explode(1) = 1;
end

h = pie(ax, counts, explode);           % h: patch,text,patch,text,...
patchH = h(1:2:end);
textH  = h(2:2:end);

% style patches (robustly in the same order as counts/labels)
for i = 1:numel(patchH)
    set(patchH(i), 'FaceColor', sliceColors(i,:), 'EdgeColor', 'w', 'LineWidth', 0.8);
end

% set text to our labels
for i = 1:numel(textH)
    set(textH(i), 'String', sliceText(i), 'Interpreter','none', ...
        'FontName', opt.FontName, 'FontSize', opt.FontSize);
end

axis(ax, 'equal');
if opt.HideAxes
    axis(ax, 'off'); % hides x/y axes
end
ax.Color = 'w';
ax.Box = 'off';

% Title
if strlength(string(opt.Title)) > 0
    title(ax, string(opt.Title), 'FontName', opt.FontName, ...
        'FontSize', opt.FontSize+2, 'FontWeight','bold');
end

% Legend
if opt.ShowLegend
    lgdStr = strings(numel(labels),1);
    for i = 1:numel(labels)
        lgdStr(i) = sprintf('%s — %.3g%% (n=%d)', labels(i), perc(i), counts(i));
    end
    legend(ax, lgdStr, 'Location','eastoutside', 'Box','off', ...
        'FontName', opt.FontName, 'FontSize', opt.FontSize);
end

% Stats output
stats = table(labels, counts, perc, 'VariableNames', {'Category','Count','Percent'});

hold(ax, 'off');
end
