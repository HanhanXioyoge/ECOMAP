function statsTbl = plot_qSub_abs_scatter(aerobicData, anaerobicData, Mul_csourcesData)
% PLOT_QSUB_ABS_SCATTER
%   Compare experimental vs model substrate uptake (q_sub) for three
%   conditions (aerobic, anaerobic, multi-carbon sources) in a single
%   scatter plot using absolute values.
%
% INPUT
%   aerobicData      : table, last two columns are [exp_q_sub, model_q_sub]
%   anaerobicData    : table, last two columns are [exp_q_sub, model_q_sub]
%   Mul_csourcesData : table, last two columns are [exp_q_sub, model_q_sub]
%
% OUTPUT
%   statsTbl : table with summary statistics per condition and for all
%              pooled data:
%              Condition, r, R2, p, N, MAE
%
% REQUIREMENT
%   This function assumes that a helper function
%       [r, p, n, mae] = compute_stats(x, y)
%   is available on the MATLAB path. R2 is computed here as r.^2.

    % ---------------------- small helper: get last 2 cols ----------------------
    getXY = @(T) deal( ...
        T{:, end-1}, ...   % x = experimental q_sub
        T{:, end});        % y = model-predicted q_sub

    % ---------------------- extract and take absolute values -------------------
    [x_a,  y_a ] = getXY(aerobicData);
    [x_an, y_an] = getXY(anaerobicData);
    [x_m,  y_m ] = getXY(Mul_csourcesData);

    % use absolute values for plotting and statistics
    x_a  = abs(x_a);   y_a  = abs(y_a);
    x_an = abs(x_an);  y_an = abs(y_an);
    x_m  = abs(x_m);   y_m  = abs(y_m);

    % ---------------------- remove NaN / Inf -----------------------------------
    valid_a  = isfinite(x_a)  & isfinite(y_a);
    valid_an = isfinite(x_an) & isfinite(y_an);
    valid_m  = isfinite(x_m)  & isfinite(y_m);

    x_a  = x_a(valid_a);   y_a  = y_a(valid_a);
    x_an = x_an(valid_an); y_an = y_an(valid_an);
    x_m  = x_m(valid_m);   y_m  = y_m(valid_m);

    % pooled data (for overall statistics)
    x_all = [x_a; x_an; x_m];
    y_all = [y_a; y_an; y_m];

    % ---------------------- statistics per condition ---------------------------
    [r_a,  p_a,  n_a,  mae_a ] = compute_stats(x_a,  y_a);
    [r_an, p_an, n_an, mae_an] = compute_stats(x_an, y_an);
    [r_m,  p_m,  n_m,  mae_m ] = compute_stats(x_m,  y_m);

    % overall statistics (pooled data)
    [r_all, p_all, n_all, mae_all] = compute_stats(x_all, y_all);

    % compute R2 as r^2
    R2_a   = r_a.^2;
    R2_an  = r_an.^2;
    R2_m   = r_m.^2;
    R2_all = r_all.^2;

    % ---------------------- build stats table ----------------------------------
    Condition = categorical( ...
        {'Aerobic'; 'Anaerobic'; 'Multi-carbon'; 'All conditions'} );

    r_vec   = [r_a;    r_an;    r_m;    r_all];
    R2_vec  = [R2_a;   R2_an;   R2_m;   R2_all];
    p_vec   = [p_a;    p_an;    p_m;    p_all];
    N_vec   = [n_a;    n_an;    n_m;    n_all];
    MAE_vec = [mae_a;  mae_an;  mae_m;  mae_all];

    statsTbl = table(Condition, r_vec, R2_vec, p_vec, N_vec, MAE_vec, ...
        'VariableNames', {'Condition','r','R2','p','N','MAE'});

    fprintf('\n==== |q_{sub}| experiment vs model: summary statistics ====\n');
    disp(statsTbl);

    % ---------------------- scatter plot (publication style) -------------------
    figure; hold on;

    % scientific color palette (Matlab default line color order)
    col_a  = [0.0000 0.4470 0.7410];  % blue   : aerobic
    col_an = [0.8500 0.3250 0.0980];  % orange : anaerobic
    col_m  = [0.4660 0.6740 0.1880];  % green  : multi-carbon

    sz = 25; % marker size

    % legend labels including r and R2
    labelA  = sprintf('Aerobic (r = %.2f, R^2 = %.2f)',  r_a,  R2_a);
    labelAn = sprintf('Anaerobic (r = %.2f, R^2 = %.2f)', r_an, R2_an);
    labelM  = sprintf('Multi C-sources (r = %.2f, R^2 = %.2f)', r_m,  R2_m);

    scatter(x_a,  y_a,  sz, 'o', ...
        'MarkerEdgeColor', col_a, ...
        'MarkerFaceColor', col_a, ...
        'MarkerFaceAlpha', 0.7, ...
        'DisplayName', labelA);

    scatter(x_an, y_an, sz, 's', ...
        'MarkerEdgeColor', col_an, ...
        'MarkerFaceColor', col_an, ...
        'MarkerFaceAlpha', 0.7, ...
        'DisplayName', labelAn);

    scatter(x_m,  y_m,  sz, '^', ...
        'MarkerEdgeColor', col_m, ...
        'MarkerFaceColor', col_m, ...
        'MarkerFaceAlpha', 0.7, ...
        'DisplayName', labelM);

    % y = x reference line in absolute value space
    allVals = [x_all; y_all];
    minVal  = min(allVals);
    maxVal  = max(allVals);

    plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 1.2, ...
        'DisplayName', 'y = x');

    xlabel('Experimental uptake rate (mmol·gDW^{-1}·h^{-1})');
    ylabel('Predicted uptake rate (mmol·gDW^{-1}·h^{-1})');

    axis equal;
    xlim([minVal maxVal]);
    ylim([minVal maxVal]);
    box on;
    grid off;  % no grid, cleaner for publication figures
    set(gcf, 'Color', 'w');

    % place legend outside to the right
    legend('Location','bestoutside');

    title('Substrate uptake (absolute values): experiment vs model');

    % annotate overall statistics on the figure (normalized axes coordinates)
    txt = sprintf('All conditions: r = %.2f, R^2 = %.2f, N = %d', ...
                  r_all, R2_all, n_all);
    text(0.02, 0.98, txt, ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontSize', 10);

end
