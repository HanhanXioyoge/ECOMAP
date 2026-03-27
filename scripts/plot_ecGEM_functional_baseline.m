% =========================================================================
% Script Name: plot_ecGEM_functional_baseline.m
% Description: Generates a grouped bar chart to assess the functional 
%              baseline of uncalibrated ecGEMs across different microbial 
%              chassis and carbon sources.
% =========================================================================

% --- Part 1: Data Input ---
% 1. Wet-lab Measurements (mu_max, h^-1)
% Order of carbon sources: Glucose, Fructose, Acetate
wet_lab_ecoli = [0.74, 0.54, 0.29];   % E. coli
wet_lab_yeast = [0.41, 0.338, 0.17];  % S. cerevisiae
wet_lab_cg    = [0.45, 0.41, 0.20];   % C. glutamicum

% 2. Uncalibrated in silico Predictions (mu_max, h^-1)
in_silico_ecoli = [0.2393, 0.2384, 0.0731]; % eciML1515
in_silico_yeast = [0.2443, 0.2511, 0.0057]; % ecYeast9
in_silico_cg    = [0.1121, 0.1123, 0.0604]; % eciCW773

% --- Part 2: Data Integration and Plotting ---
% Concatenate block data into full column vectors (9x1)
wet_lab = [wet_lab_ecoli, wet_lab_yeast, wet_lab_cg]';
in_silico = [in_silico_ecoli, in_silico_yeast, in_silico_cg]';

% Combine into a matrix for the grouped bar chart
Y = [wet_lab, in_silico];

% Create the figure window
figure('Name', 'ecGEM Functional Baseline', 'Position', [100, 100, 1100, 600], 'Color', 'w');

% Plot grouped bars
b = bar(Y, 'grouped');

% Set bar colors: Dark grey-blue (Wet-lab), Light blue (In silico)
b(1).FaceColor = [0.35, 0.40, 0.45]; 
b(1).EdgeColor = 'k';
b(2).FaceColor = [0.60, 0.75, 0.90]; 
b(2).EdgeColor = 'k';

hold on;

% Customize axes and labels
set(gca, 'XTick', 1:9, 'FontSize', 11);
xticklabels({'Glucose', 'Fructose', 'Acetate', ...
             'Glucose', 'Fructose', 'Acetate', ...
             'Glucose', 'Fructose', 'Acetate'});
         
ylabel('Maximum Specific Growth Rate, \mu_{max} (h^{-1})', 'FontSize', 13, 'FontWeight', 'bold');
xlabel('Microbial Chassis and Carbon Source', 'FontSize', 13, 'FontWeight', 'bold');

% Set Y-axis limits dynamically based on the maximum wet-lab value
max_y = max(wet_lab) * 1.3;
ylim([0, max_y]); 

% Add legend
legend({'Wet-lab Measurement', 'Uncalibrated in silico Prediction'}, ...
    'Location', 'northwest', 'FontSize', 11);

% Add chassis separation lines
xline(3.5, '--k', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
xline(6.5, '--k', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);

% Add chassis text labels dynamically based on max height
text_height = max_y * 0.92;
text(2, text_height, 'E. coli', 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontAngle', 'italic', 'FontWeight', 'bold');
text(5, text_height, 'S. cerevisiae', 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontAngle', 'italic', 'FontWeight', 'bold');
text(8, text_height, 'C. glutamicum', 'HorizontalAlignment', 'center', 'FontSize', 15, 'FontAngle', 'italic', 'FontWeight', 'bold');

% Add percentage labels (fractional coverage) above the in silico bars
[ngroups, nbars] = size(Y);
x_coords = nan(ngroups, nbars);
for i = 1:nbars
    x_coords(:,i) = b(i).XEndPoints;
end

percentages = (in_silico ./ wet_lab) * 100;
for i = 1:9
    text(x_coords(i,2), Y(i,2) + (max_y * 0.03), sprintf('~%.1f%%', percentages(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k');
end

% Final touches: title and grid
title('Functional Baseline Assessment of Uncalibrated ecGEMs', 'FontSize', 15, 'FontWeight', 'bold');
grid on;
set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.6);
hold off;