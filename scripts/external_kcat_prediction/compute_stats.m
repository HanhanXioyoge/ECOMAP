function [r, p, n, mae] = compute_stats(x, y)
% Compute Pearson r, p-value, number of points, and MAE in log10 space.
    if nargin < 2, r = NaN; p = NaN; n = 0; mae = NaN; return; end
    if isempty(x) || isempty(y), r = NaN; p = NaN; n = 0; mae = NaN; return; end
    % keep finite
    finiteMask = isfinite(x) & isfinite(y);
    x = x(finiteMask); y = y(finiteMask);
    n = numel(x);
    if n >= 2
        C = corrcoef(x, y, 'Rows','pairwise');
        r = C(1,2);
        t = r * sqrt((n-2) / max(1e-12, 1 - r^2));
        try
            p = 2 * tcdf(-abs(t), n-2);
        catch
            p = NaN;
        end
        mae = mean(abs(x - y));
    elseif n == 1
        r = NaN; p = NaN; mae = abs(x - y);
    else
        r = NaN; p = NaN; mae = NaN;
    end
end