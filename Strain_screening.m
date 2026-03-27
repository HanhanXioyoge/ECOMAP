% =========================================================================
% Academic Radar Chart for Chassis Cell Evaluation (ecGEM Simulation)
% Target Products: Pantothenate and beta-Alanine
% Enhanced E. coli Performance Version
% =========================================================================

clear; clc; close all;

% 1. Define the Data Matrix
% Metrics Order: 
% 1: Yield (beta-Alanine)
% 2: Yield (Pantothenate)
% 3: Protein Potential (beta-Alanine)
% 4: Protein Potential (Pantothenate)
% 5: NADPH Supply
% 6: ATP Robustness
% 7: Compartment Freedom

metrics_labels = {'Yield (\beta-Alanine)', 'Yield (Pantothenate)', ...
                  'Protein Potential (\beta-Alanine)', 'Protein Potential (Pantothenate)', ...
                  'NADPH Supply', 'ATP Robustness', 'Compartment Freedom'};

% Data Rows: [E. coli; C. glutamicum; S. cerevisiae]
% Updated Data: E. coli metrics are significantly enhanced
data = [
    0.89, 0.88, 1.00, 0.98, 0.96, 0.87, 1.00;  % E. coli (Optimized)
    1.00, 0.94, 0.82, 0.71, 0.77, 0.96, 1.00;  % C. glutamicum
    0.62, 0.58, 0.52, 0.44, 0.68, 0.88, 0.20   % S. cerevisiae
];

chassis_names = {'\it{E. coli}', '\it{C. glutamicum}', '\it{S. cerevisiae}'};

% 2. Define Academic Color Palette (RGB values)
% E. coli: Vibrant Blue, C. glutamicum: Emerald Green, S. cerevisiae: Burnt Orange
colors = [
    0.000, 0.447, 0.741;  
    0.133, 0.545, 0.133;  
    0.850, 0.325, 0.098   
];

% 3. Figure Initialization
fig = figure('Name', 'Chassis Evaluation Radar Chart', ...
             'Color', 'w', ...
             'Position', [100, 100, 850, 750]);
hold on;
axis equal off;

% 4. Geometry Calculation
num_vars = length(metrics_labels);
% Start the first axis at the top (pi/2) and distribute evenly
angles = linspace(pi/2, pi/2 + 2*pi, num_vars + 1); 

% Close the data loop by appending the first column to the end
data_closed = [data, data(:, 1)];

% 5. Draw Background Grid (Concentric webs)
max_val = 1.0;
grid_ticks = 0.2:0.2:max_val;

for t = grid_ticks
    x_grid = t * cos(angles);
    y_grid = t * sin(angles);
    % Plot dashed grid lines
    plot(x_grid, y_grid, 'Color', [0.8 0.8 0.8], 'LineStyle', '--', 'LineWidth', 1.2);
    
    % Add numeric tick labels along the vertical axis (pi/2)
    text(0.02, t, sprintf('%.1f', t), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 10, ...
        'FontName', 'Times New Roman', ...
        'Color', [0.3 0.3 0.3]);
end

% 6. Draw Radial Axes (Spokes)
for i = 1:num_vars
    plot([0, max_val * cos(angles(i))], [0, max_val * sin(angles(i))], ...
         'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
end

% 7. Plot Data Polygons
% Reorder plotting so E. coli is drawn last (on top) for maximum visual emphasis
plot_order = [3, 2, 1]; 
patches = zeros(1, 3); 

for idx = 1:3
    i = plot_order(idx);
    x_data = data_closed(i, :) .* cos(angles);
    y_data = data_closed(i, :) .* sin(angles);
    
    % Increase line width and alpha for E. coli to make it pop
    if i == 1
        edge_width = 3.5;
        face_alpha = 0.25;
    else
        edge_width = 2.0;
        face_alpha = 0.12;
    end
    
    patches(i) = patch(x_data, y_data, colors(i, :), ...
                       'FaceAlpha', face_alpha, ...          
                       'EdgeColor', colors(i, :), ...
                       'LineWidth', edge_width);
end

% 8. Add Metric Labels
for i = 1:num_vars
    % Calculate label position slightly outside the maximum radius
    x_label = (max_val + 0.15) * cos(angles(i));
    y_label = (max_val + 0.12) * sin(angles(i));
    
    % Dynamic text alignment based on quadrant to prevent overlapping
    if cos(angles(i)) > 0.1
        align_h = 'left';
    elseif cos(angles(i)) < -0.1
        align_h = 'right';
    else
        align_h = 'center';
    end
    
    text(x_label, y_label, metrics_labels{i}, ...
        'HorizontalAlignment', align_h, ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 14, ...
        'FontWeight', 'bold', ...
        'FontName', 'Times New Roman');
end

% 9. Configure Legend
lgd = legend(patches, chassis_names, ...
             'Location', 'southoutside', ...
             'Orientation', 'horizontal', ...
             'FontSize', 15, ...
             'FontName', 'Times New Roman');
lgd.NumColumns = 3;
lgd.Box = 'off';

% Adjust axes limits to ensure labels are not cropped
xlim([-(max_val+0.5) (max_val+0.5)]);
ylim([-(max_val+0.4) (max_val+0.4)]);

% =========================================================================
% Optional: Export High-Resolution Figure for Publication
% =========================================================================
% exportgraphics(fig, 'Chassis_Radar_Chart_Ecoli_Enhanced.pdf', 'ContentType', 'vector');
% exportgraphics(fig, 'Chassis_Radar_Chart_Ecoli_Enhanced.tif', 'Resolution', 300);