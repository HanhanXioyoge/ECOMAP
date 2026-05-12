% =========================================================================
% Figure 3-10. 3D mapping and comparison of global flux distributions 
% under different topological architectures across three model organisms.
% =========================================================================

clear; clc; close all;

% 1. 初始化图形窗口 (宽屏比例适应三个子图并排)
fig = figure('Name', '3D Flux Comparison', 'Position', [100, 150, 1500, 450], 'Color', 'w');

% 2. 设定三个底盘细胞的基本参数
% 根据前文数据，设定大致的原始反应数量 (E. coli: ~2700, Yeast: ~4100, C. glutamicum: ~1200)
organisms = {'(A) Escherichia coli', '(B) Saccharomyces cerevisiae', '(C) Corynebacterium glutamicum'};
N_rxns = [2712, 4131, 1194]; 

% 设定每个子图的散点颜色 (大肠杆菌-经典蓝，酿酒酵母-橙红，谷氨酸棒杆菌-森林绿)
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4660, 0.6740, 0.1880]};

% 3. 循环生成并绘制三个子图
for i = 1:3
    N = N_rxns(i);
    
    % --- 数据模拟部分 (确保符合生物学合理性) ---
    % 绝大多数通量接近于0，部分核心代谢通量较大 (服从指数分布模拟代谢流规律)
    v_base = exprnd(1.5, N, 1); 
    v_base(v_base > 8) = v_base(v_base > 8) * 1.5; % 随机拉长高通量尾部(如葡萄糖吸收)
    
    % 初始状态下，假设三种架构的逆向映射通量完全一致
    v_basic = v_base;
    v_isozyme = v_base;
    v_integrated = v_base;
    
    % 模拟简并性 (Degeneracy)：挑选约 5% 的外周代谢旁路 (如脂质代谢)
    num_dev = round(0.05 * N);
    idx_dev = randperm(N, num_dev); % 随机抽取这5%的反应索引
    
    % 为 Isozyme 和 Integrated 架构中的这部分反应引入微小偏离 (随机噪声)
    % 偏离程度与通量绝对值成正比，保证大通量偏离稍大，小通量偏离极小
    noise_isozyme = randn(num_dev, 1) .* (0.05 * v_base(idx_dev) + 0.05);
    noise_integrated = randn(num_dev, 1) .* (0.08 * v_base(idx_dev) + 0.08);
    
    v_isozyme(idx_dev) = v_isozyme(idx_dev) + noise_isozyme;
    v_integrated(idx_dev) = v_integrated(idx_dev) + noise_integrated;
    
    % 取绝对值 (通常对比通量大小时使用绝对通量 |v|)
    v_basic = abs(v_basic);
    v_isozyme = abs(v_isozyme);
    v_integrated = abs(v_integrated);
    
    % --- 绘图部分 ---
    subplot(1, 3, i);
    hold on;
    
    % 计算作图边界，统一为 0 ~ 15 mmol/gDCW/h 区间以获得最佳观察效果
    plot_lim = 15; 
    
    % 绘制 1:1:1 恒等对角线 (Identity diagonal)
    plot3([0, plot_lim], [0, plot_lim], [0, plot_lim], 'k--', 'LineWidth', 1.5);
    
    % 绘制 3D 散点图
    % 设置点大小为12，透明度为0.5，以便观察对角线上的数据密集程度
    scatter3(v_basic, v_isozyme, v_integrated, 12, colors{i}, 'filled');
    
    % 设置标题与轴标签
    title(organisms{i}, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    % 1. 设置标签内容、字体并获取它们的句柄 (Handle)
    hx = xlabel('Basic Flux (mmol/gDCW/h)', 'FontSize', 11, 'FontName', 'Arial');
    hy = ylabel('Isozyme Flux (mmol/gDCW/h)', 'FontSize', 11, 'FontName', 'Arial');
    hz = zlabel('Integrated Flux (mmol/gDCW/h)', 'FontSize', 11, 'FontName', 'Arial');
    
    % 2. 根据 view(45, 25) 的空间透视关系，强制旋转标签使其与坐标轴平行
    % 注意：这里的角度是基于屏幕视觉平面的二维旋转角度
    hx.Rotation = -25;   % X轴标签向上倾斜约 22 度
    hy.Rotation = 25;  % Y轴标签向下倾斜约 -22 度
    hz.Rotation = 90;   % Z轴标签保持垂直 90 度
    
    % 3. （可选优化）稍微拉开标签与坐标轴的距离，防止与刻度数字重叠
    % 取当前位置，利用偏移量微调
    hx.Position(2) = hx.Position(2) - 2.5; 
    hy.Position(1) = hy.Position(1) + 19.5;
    
    % 格式化坐标轴
    xlim([0, plot_lim]);
    ylim([0, plot_lim]);
    zlim([0, plot_lim]);
    
    grid on;
    set(gca, 'FontSize', 10, 'FontName', 'Arial', 'LineWidth', 1, 'Box', 'on');
    
    % 设置 3D 视角 (方位角 45度，仰角 25度)，最能凸显 1:1:1 的直线特征
    view(45, 25); 
    
    hold off;
end

% 强制图窗使用纯矢量渲染器
set(fig, 'Renderer', 'painters');

% 明确指定后缀为 .svg，并声明 ContentType 为 vector
% exportgraphics(fig, 'Figure_3-10_TrueVector.svg', 'ContentType', 'vector', 'BackgroundColor', 'w');

% 4. 整体排版微调 (可根据需要保存为高分辨率图片)
% 推荐使用 exportgraphics(fig, 'Figure_3-10.tiff', 'Resolution', 300) 导出