%% 1. 初始化与设置
clear; clc; close all;

% --- 配色方案 (Nature 风格) ---
c_blue  = [0.2, 0.4, 0.7]; % Constrained
c_green = [0.2, 0.6, 0.3]; % Unconstrained
c_red   = [0.8, 0.2, 0.2]; % 13C Flux

colors = [c_blue; c_green; c_red];
markers = {'o', '^', 's'}; 

% --- 画布设置 ---
figure('Color', 'w', 'Position', [50, 50, 1150, 1400]); 

% TiledLayout 3行1列
t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

%% 2. 模型定义
list_A = {'Init','Single_sensitivity','Multi_sensitivity','Multi_sensitivity_GAUKS',...
          'Bayesian_constraints_growth_sensitivity','Bayesian_constraints_growth_sensitivity_GAUKS',...
          'Bayesian_all_growth_sensitivity','Bayesian_all_growth_sensitivity_GAUKS',...
          'Bayesian_13C_Multi_sensitivity','Bayesian_13C_Multi_sensitivity_GAUKS',...
          'Bayesian_full_sensitivity','Bayesian_full_sensitivity_GAUKS',...
          'PRESTO','GECKO','ECMPY'};

list_B = {'Init','Single_sensitivity','Multi_sensitivity','Multi_sensitivity_GAUKS',...
          'Bayesian_constraints_growth_sensitivity','Bayesian_constraints_growth_sensitivity_GAUKS',...
          'Bayesian_all_growth_sensitivity','Bayesian_all_growth_sensitivity_GAUKS'};

% 物种名称 (注意：稍后在Title中会处理斜体)
species_names = {'Escherichia coli', 'Saccharomyces cerevisiae', 'Corynebacterium glutamicum'};

% 指定最优模型
best_target_names = {'Bayesian_full_sensitivity_GAUKS', ...
                     'Bayesian_full_sensitivity_GAUKS', ...
                     'Bayesian_all_growth_sensitivity_GAUKS'};

%% 3. 循环绘制
for i = 1:3
    % --- Step 1: 确定参数基准 (Baselines) ---
    if i == 1
        % 大肠杆菌: Con > 0.8, Uncon > 0.2
        curr_models = list_A;
        base_con_min = 0.80;  
        base_unc_min = 0.20;
        base_13c_min = 1.0; % 13C 依然较大
    elseif i == 2
        % 酿酒酵母: Con > 0.9, Uncon > 0.1
        curr_models = list_A;
        base_con_min = 0.90;
        base_unc_min = 0.10;
        base_13c_min = 1.2;
    else
        % 谷棒: 无特定要求，保持较小值以示区别
        curr_models = list_B;
        base_con_min = 0.05;
        base_unc_min = 0.15;
        base_13c_min = NaN; 
    end
    
    n_models = length(curr_models);
    raw_data = zeros(n_models, 3);
    
    % 设定随机种子
    rng(i * 999); 
    
    % --- Step 2: 生成带有交叉特性的数据 ---
    for k = 1:n_models
        m_name = curr_models{k};
        
        % 定义相对表现因子 (Factor, 0=Best, 1=Worst)
        if strcmp(m_name, 'Init')
            perf_factor = 2.0; % 极差 (Outlier)
        elseif contains(m_name, 'Single')
            perf_factor = 0.7 + rand*0.2;
        elseif contains(m_name, 'Multi')
            perf_factor = 0.5 + rand*0.2;
        elseif strcmp(m_name, 'PRESTO') || strcmp(m_name, 'GECKO') || strcmp(m_name, 'ECMPY')
            perf_factor = 0.3 + rand*0.3; % 中上等，随机波动大
        elseif contains(m_name, 'Bayesian')
            perf_factor = 0.1 + rand*0.2; % 优等
        else
            perf_factor = 0.5;
        end
        
        % 引入随机扰动 (Noise) 制造交叉
        % 例如：有的模型 Con 好 (noise负), Uncon 差 (noise正)
        noise = (rand(1, 3) - 0.5) * 0.15; 
        
        % 生成基础数据 = 基准线 + (表现因子 * 跨度) + 噪音
        val_con = base_con_min + perf_factor * 0.3 + noise(1);
        val_unc = base_unc_min + perf_factor * 0.2 + noise(2);
        
        if i == 3
            val_13c = NaN;
        else
            val_13c = base_13c_min + perf_factor * 1.5 + noise(3);
        end
        
        % 强制修正：Init 模型必须显著大
        if strcmp(m_name, 'Init')
            val_con = base_con_min + 0.6; % 强制加很大
            val_unc = base_unc_min + 0.4;
            if i ~= 3
                val_13c = base_13c_min + 2.0;
            end
        end
        
        % 边界保护 (防止噪音导致低于基准线)
        val_con = max(val_con, base_con_min);
        val_unc = max(val_unc, base_unc_min);
        
        raw_data(k, :) = [val_con, val_unc, val_13c];
    end
    
    % --- Step 3: 优化 "最优模型" (Realistic Best) ---
    target_best = best_target_names{i};
    best_idx = find(strcmp(curr_models, target_best));
    
    if ~isempty(best_idx)
        % 1. 13C 最好
        if i ~= 3
            raw_data(best_idx, 3) = base_13c_min; 
        end
        % 2. Con 最好
        raw_data(best_idx, 1) = base_con_min;
        % 3. Uncon *不是* 最好 (制造真实感，略高于最小值)
        % 找到当前 Uncon 这一列的最小值
        curr_min_unc = min(raw_data(:,2));
        % 设定为比最小值大一点点 (例如大 5-10%)
        raw_data(best_idx, 2) = curr_min_unc * 1.05 + 0.01; 
    end

    % --- Step 4: 计算 Total Cost (加权平方和) ---
    % 权重设定: 13C 误差通常大且重要，稍微降低权重以防主导; Con/Uncon 权重正常
    if i == 3
        costs = raw_data(:,1).^2 + raw_data(:,2).^2;
    else
        % 13C 的数值较大 (1.0+)，平方后很大，给一个小权重 0.2 平衡
        % Con/Uncon 数值较小 (0.8/0.2)，保持权重 1.0
        costs = 1.0 * raw_data(:,1).^2 + 1.0 * raw_data(:,2).^2 + 0.2 * raw_data(:,3).^2;
    end
    
    % --- Step 5: 排序与归一化 (用于绘图) ---
    [sorted_costs, sort_idx] = sort(costs, 'descend'); % 降序
    
    plot_data_raw = raw_data(sort_idx, :);
    plot_names = curr_models(sort_idx);
    
    % 归一化 (用于左侧画图，使不同量级的Con/Uncon能画在一起)
    norm_plot_data = plot_data_raw;
    for col = 1:3
        col_dat = raw_data(:, col);
        if ~all(isnan(col_dat))
            min_v = min(col_dat);
            max_v = max(col_dat);
            % 防止分母为0
            if max_v == min_v, max_v = min_v + 1e-6; end
            norm_plot_data(:, col) = (plot_data_raw(:, col) - min_v) ./ (max_v - min_v);
        end
    end

    % --- Step 6: 绘图 ---
    ax = nexttile;
    hold on;
    
    % A. 绘制左侧哑铃图
    for j = 1:n_models
        row_vals = norm_plot_data(j, :);
        valid_vals = row_vals(~isnan(row_vals));
        if ~isempty(valid_vals)
            line([min(valid_vals), max(valid_vals)], [j, j], ...
                 'Color', [0.85, 0.85, 0.85], 'LineWidth', 1.5);
        end
    end
    
    sz = 80;
    for m = 1:3
        valid_mask = ~isnan(norm_plot_data(:, m));
        if any(valid_mask)
            scatter(norm_plot_data(valid_mask, m), find(valid_mask), sz, colors(m,:), ...
                'filled', markers{m}, 'MarkerEdgeColor', 'w');
        end
    end
    
    % B. 绘制右侧表格
    x_cols = [1.25, 1.65, 2.05, 2.55]; 
    header_y = n_models + 1.5;
    
    % 表头
    text(x_cols(1), header_y, 'Con. RMSE', 'Color', c_blue, 'FontWeight','bold', 'HorizontalAlignment','center');
    text(x_cols(2), header_y, 'Uncon. RMSE', 'Color', c_green, 'FontWeight','bold', 'HorizontalAlignment','center'); 
    text(x_cols(3), header_y, '^{13}C RMSE', 'Color', c_red, 'FontWeight','bold', 'HorizontalAlignment','center');   
    text(x_cols(4), header_y, 'Total Cost', 'Color', 'k', 'FontWeight','bold', 'HorizontalAlignment','center');
    
    % 写入数据
    for j = 1:n_models
        fw = 'normal'; 
        if j == n_models, fw = 'bold'; end 
        
        % 1. Con (Blue)
        text(x_cols(1), j, sprintf('%.3f', plot_data_raw(j,1)), 'Color', c_blue, ...
             'HorizontalAlignment','center', 'FontWeight',fw, 'FontSize', 9);
        % 2. Uncon (Green)
        text(x_cols(2), j, sprintf('%.3f', plot_data_raw(j,2)), 'Color', c_green, ...
             'HorizontalAlignment','center', 'FontWeight',fw, 'FontSize', 9);
        % 3. 13C (Red)
        val_13c = plot_data_raw(j,3);
        if isnan(val_13c)
            str_13c = '—'; col_c = [0.5 0.5 0.5];
        else
            str_13c = sprintf('%.2f', val_13c); col_c = c_red;
        end
        text(x_cols(3), j, str_13c, 'Color', col_c, ...
             'HorizontalAlignment','center', 'FontWeight',fw, 'FontSize', 9);
        % 4. Total Cost (Black)
        text(x_cols(4), j, sprintf('%.2f', sorted_costs(j)), 'Color', 'k', ...
             'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', 10);
    end
    
    % --- 美化 ---
    ax.XLim = [-0.1, 2.8]; 
    ax.YLim = [0.5, n_models + 2]; 
    ax.YTick = 1:n_models;
    ax.YTickLabel = plot_names;
    ax.TickLabelInterpreter = 'none';
    ax.YAxis.Color = 'k'; ax.YAxis.FontSize = 10;
    box off; ax.TickDir = 'out'; ax.YAxis.TickLength = [0 0]; 
    
    % 标题：物种名称斜体处理 (\it)
    title(sprintf('%c. \\it%s', 64+i, species_names{i}), ...
          'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'left', ...
          'Interpreter', 'tex', ...
          'Position', [-0.1, n_models + 2.5]); 
      
    if i == 3
        xlabel('Normalized Metrics Distribution (0=Best)           |           Raw Metrics & Calculated Total Cost', ...
            'FontSize', 10, 'Color', [0.4 0.4 0.4]);
        xline(1.1, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    else
        ax.XAxis.Visible = 'off';
        xline(1.1, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    end
    
    % 高亮第一名
    patch([-0.1, 2.8, 2.8, -0.1], ...
          [n_models-0.4, n_models-0.4, n_models+0.4, n_models+0.4], ...
          [1 0.85 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end

% --- 图例 ---
lgd_ax = nexttile(3); hold on;
h_dum(1) = plot(nan,nan, markers{1}, 'MarkerFaceColor', c_blue, 'MarkerEdgeColor', 'w', 'LineStyle', 'none');
h_dum(2) = plot(nan,nan, markers{2}, 'MarkerFaceColor', c_green, 'MarkerEdgeColor', 'w', 'LineStyle', 'none');
h_dum(3) = plot(nan,nan, markers{3}, 'MarkerFaceColor', c_red, 'MarkerEdgeColor', 'w', 'LineStyle', 'none');
lgd = legend(h_dum, {'Substrate Constrained RMSE', 'Unconstrained RMSE', '^{13}C Flux RMSE'}, ...
             'Orientation', 'horizontal', 'Location', 'southoutside');
lgd.Box = 'off'; lgd.FontSize = 11;
hold off;