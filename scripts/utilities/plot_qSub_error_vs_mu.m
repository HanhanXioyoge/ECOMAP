function statsTbl = plot_qSub_error_vs_mu(aerobicData, anaerobicData, Mul_csourcesData)
% PLOT_QSUB_ERROR_VS_MU
%   Plot prediction error of substrate uptake vs the first-column variable
%   (e.g. growth rate mu) for three condition sets:
%   aerobic, anaerobic, and multi-carbon sources.
%
%   error = |q_model| - |q_exp|
%
%   NOTE: The reported r, R2, MAE are computed between |q_exp| and
%   |q_model| (goodness-of-fit of the model), but they are annotated
%   on the error-vs-mu plot for convenience.
%
% INPUT
%   aerobicData      : table, first column = Z variable (e.g. mu),
%                      last two columns = [exp_q_sub, model_q_sub]
%   anaerobicData    : same structure
%   Mul_csourcesData : same structure
%
% OUTPUT
%   statsTbl : table with summary statistics between |q_exp| and |q_model|
%              per condition and for all pooled data:
%              Condition, r, R2, p, N, MAE

    % ---------------------- helper: get X, Y, Z -------------------------------
    % X = experimental q_sub (next-to-last column)
    % Y = model-predicted q_sub (last column)
    % Z = first-column variable (e.g. growth rate mu)
    getXYZ = @(T) deal( ...
        T{:, end-1}, ...   % x = experimental q_sub
        T{:, end},   ...   % y = model-predicted q_sub
        T{:, 1});          % z = first-column variable (e.g. mu)

    % ---------------------- extract raw data ----------------------------------
    [x_a,  y_a,  z_a ] = getXYZ(aerobicData);
    [x_an, y_an, z_an] = getXYZ(anaerobicData);
    [x_m,  y_m,  z_m ] = getXYZ(Mul_csourcesData);

    % use absolute values for q_sub when comparing model vs experiment
    x_a  = abs(x_a);   y_a  = abs(y_a);
    x_an = abs(x_an);  y_an = abs(y_an);
    x_m  = abs(x_m);   y_m  = abs(y_m);

    % ---------------------- remove NaN / Inf ----------------------------------
    valid_a  = isfinite(x_a)  & isfinite(y_a)  & isfinite(z_a);
    valid_an = isfinite(x_an) & isfinite(y_an) & isfinite(z_an);
    valid_m  = isfinite(x_m)  & isfinite(y_m)  & isfinite(z_m);

    x_a  = x_a(valid_a);   y_a  = y_a(valid_a);   z_a  = z_a(valid_a);
    x_an = x_an(valid_an); y_an = y_an(valid_an); z_an = z_an(valid_an);
    x_m  = x_m(valid_m);   y_m  = y_m(valid_m);   z_m  = z_m(valid_m);

    % ---------------------- prediction error ----------------------------------
    err_a  = y_a  - x_a;
    err_an = y_an - x_an;
    err_m  = y_m  - x_m;

    % pooled data (for overall fit between |q_exp| and |q_model|)
    x_all   = [x_a;  x_an;  x_m];
    y_all   = [y_a;  y_an;  y_m];
    err_all = [err_a; err_an; err_m];

    % ---------------------- stats: |q_exp| vs |q_model| -----------------------
    [r_a,  p_a,  n_a,  mae_a ] = compute_stats(x_a,  y_a);
    [r_an, p_an, n_an, mae_an] = compute_stats(x_an, y_an);
    [r_m,  p_m,  n_m,  mae_m ] = compute_stats(x_m,  y_m);
    [r_all, p_all, n_all, mae_all] = compute_stats(x_all, y_all);

    R2_a   = r_a.^2;
    R2_an  = r_an.^2;
    R2_m   = r_m.^2;
    R2_all = r_all.^2;

    % ---------------------- build and display stats table ---------------------
    Condition = categorical( ...
        {'Aerobic'; 'Anaerobic'; 'Multi-carbon'; 'All conditions'} );

    r_vec   = [r_a;    r_an;    r_m;    r_all];
    R2_vec  = [R2_a;   R2_an;   R2_m;   R2_all];
    p_vec   = [p_a;    p_an;    p_m;    p_all];
    N_vec   = [n_a;    n_an;    n_m;    n_all];
    MAE_vec = [mae_a;  mae_an;  mae_m;  mae_all];

    statsTbl = table(Condition, r_vec, R2_vec, p_vec, N_vec, MAE_vec, ...
        'VariableNames', {'Condition','r','R2','p','N','MAE'});

    fprintf('\n==== |q_{sub}| experiment vs model: fit statistics ====\n');
    disp(statsTbl);

    % ---------------------- plotting: error vs mu -----------------------------
    figure('Color','w'); hold on;

    % colors
    col_a  = [0.0000 0.4470 0.7410];  % blue   : aerobic
    col_an = [0.8500 0.3250 0.0980];  % orange : anaerobic
    col_m  = [0.4660 0.6740 0.1880];  % green  : multi-carbon

    sz = 25; % marker size

    scatter(z_a,  err_a,  sz, 'o', ...
        'MarkerEdgeColor', col_a, ...
        'MarkerFaceColor', col_a, ...
        'MarkerFaceAlpha', 0.7, ...
        'DisplayName', 'Aerobic');

    scatter(z_an, err_an, sz, 's', ...
        'MarkerEdgeColor', col_an, ...
        'MarkerFaceColor', col_an, ...
        'MarkerFaceAlpha', 0.7, ...
        'DisplayName', 'Anaerobic');

    scatter(z_m,  err_m,  sz, '^', ...
        'MarkerEdgeColor', col_m, ...
        'MarkerFaceColor', col_m, ...
        'MarkerFaceAlpha', 0.7, ...
        'DisplayName', 'Multi C-sources');

    % zero-error reference line
    yline(0, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Zero error');

    xlabel('Growth rate \mu (h^{-1})');  % adjust if Z is not mu
    ylabel('Prediction error (|q_{model}| - |q_{exp}|) (mmol·gDW^{-1}·h^{-1})');

    box on;
    legend('Location','best');

    title('Prediction error in substrate uptake vs growth rate');

    % annotate fit statistics (between |q_exp| and |q_model|) on the figure
    txt = sprintf([ ...
        'Aerobic:      r = %.2f, R^2 = %.2f, N = %d, MAE = %.2f\n' ...
        'Anaerobic:    r = %.2f, R^2 = %.2f, N = %d, MAE = %.2f\n' ...
        'Multi C-srcs: r = %.2f, R^2 = %.2f, N = %d, MAE = %.2f\n' ...
        'All cond.:    r = %.2f, R^2 = %.2f, N = %d, MAE = %.2f'], ...
        r_a,   R2_a,   n_a,   mae_a, ...
        r_an,  R2_an,  n_an,  mae_an, ...
        r_m,   R2_m,   n_m,   mae_m, ...
        r_all, R2_all, n_all, mae_all);

    text(0.02, 0.98, txt, ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontSize', 9);

end
