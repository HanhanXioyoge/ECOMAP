% =========================================================================
%
% Script Name: plot_ecGEM_functional_baseline.m
%
% Description: Generates a grouped bar chart to assess the functional
%              baseline and calibration performance of ecGEMs across
%              different microbial chassis and carbon sources.
%
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


% 3. Calibrated in silico Predictions (mu_max, h^-1)
%
% NOTE:
% The following values are PLACEHOLDER values for visualization only.
% TODO: Replace them with the actual calibrated growth rates.
%
% Order of carbon sources: Glucose, Fructose, Acetate

calibrated_ecoli = [0.68, 0.50, 0.25];   % PLACEHOLDER - TODO: replace
calibrated_yeast = [0.38, 0.31, 0.14];   % PLACEHOLDER - TODO: replace
calibrated_cg    = [0.41, 0.37, 0.17];   % PLACEHOLDER - TODO: replace


% --- Part 2: Data Integration and Plotting ---

% Concatenate block data into full column vectors (9x1)

wet_lab = [wet_lab_ecoli, wet_lab_yeast, wet_lab_cg]';

in_silico = [in_silico_ecoli, in_silico_yeast, in_silico_cg]';

calibrated = [calibrated_ecoli, calibrated_yeast, calibrated_cg]';


% Combine into a matrix for the grouped bar chart
% Column 1: Wet-lab
% Column 2: Uncalibrated prediction
% Column 3: Calibrated prediction

Y = [wet_lab, in_silico, calibrated];


% Create the figure window

figure('Name', 'ecGEM Functional Baseline', ...
       'Position', [100, 100, 1150, 620], ...
       'Color', 'w');


% Plot grouped bars

b = bar(Y, 'grouped');

hold on;


% Set bar colors
%
% Wet-lab:       #4f5c6b
% Uncalibrated:  #91aecf
% Calibrated:    #6e535e

color_wet_lab      = [79, 92, 107] / 255;   % #4f5c6b
color_uncalibrated = [145, 174, 207] / 255; % #91aecf
color_calibrated   = [110, 83, 94] / 255;   % #6e535e

b(1).FaceColor = color_wet_lab;
b(1).EdgeColor = [0.20, 0.20, 0.20];
b(1).LineWidth = 0.8;

b(2).FaceColor = color_uncalibrated;
b(2).EdgeColor = [0.20, 0.20, 0.20];
b(2).LineWidth = 0.8;

b(3).FaceColor = color_calibrated;
b(3).EdgeColor = [0.20, 0.20, 0.20];
b(3).LineWidth = 0.8;


% Customize axes and labels

set(gca, ...
    'XTick', 1:9, ...
    'FontSize', 11, ...
    'LineWidth', 1.0, ...
    'Box', 'on');

xticklabels({'Glucose', 'Fructose', 'Acetate', ...
             'Glucose', 'Fructose', 'Acetate', ...
             'Glucose', 'Fructose', 'Acetate'});

ylabel('Maximum Specific Growth Rate, \mu_{max} (h^{-1})', ...
       'FontSize', 13, ...
       'FontWeight', 'bold');

xlabel('Microbial Chassis and Carbon Source', ...
       'FontSize', 13, ...
       'FontWeight', 'bold');


% Set Y-axis limits dynamically based on the maximum value

max_y_value = max(Y(:));
max_y = max_y_value * 1.30;

ylim([0, max_y]);


% Add legend with white background and visible border

lgd = legend({'Wet-lab Measurement', ...
              'Uncalibrated in silico Prediction', ...
              'Calibrated in silico Prediction'}, ...
             'Location', 'northwest', ...
             'FontSize', 10.5);

lgd.Box = 'on';
lgd.Color = 'w';
lgd.EdgeColor = [0.35, 0.35, 0.35];
lgd.LineWidth = 0.8;


% Add chassis separation lines

xline(3.5, '--', ...
      'LineWidth', 1.3, ...
      'Color', [0.55, 0.55, 0.55]);

xline(6.5, '--', ...
      'LineWidth', 1.3, ...
      'Color', [0.55, 0.55, 0.55]);


% Add chassis text labels dynamically based on max height

text_height = max_y * 0.93;

text(2, text_height, 'E. coli', ...
     'HorizontalAlignment', 'center', ...
     'FontSize', 14, ...
     'FontAngle', 'italic', ...
     'FontWeight', 'bold');

text(5, text_height, 'S. cerevisiae', ...
     'HorizontalAlignment', 'center', ...
     'FontSize', 14, ...
     'FontAngle', 'italic', ...
     'FontWeight', 'bold');

text(8, text_height, 'C. glutamicum', ...
     'HorizontalAlignment', 'center', ...
     'FontSize', 14, ...
     'FontAngle', 'italic', ...
     'FontWeight', 'bold');


% Obtain exact x-coordinates of individual bars

[ngroups, nbars] = size(Y);

x_coords = nan(ngroups, nbars);

for i = 1:nbars
    x_coords(:,i) = b(i).XEndPoints;
end


% Calculate percentage of wet-lab growth rate
%
% Uncalibrated coverage:
% predicted growth rate / wet-lab growth rate

uncalibrated_percentages = (in_silico ./ wet_lab) * 100;


% Calibrated coverage:
% calibrated growth rate / wet-lab growth rate
%
% NOTE:
% These percentages currently depend on the PLACEHOLDER calibrated values.

calibrated_percentages = (calibrated ./ wet_lab) * 100;


% Add percentage labels above the uncalibrated bars

for i = 1:9

    text(x_coords(i,2), ...
         Y(i,2) + (max_y * 0.020), ...
         sprintf('%.1f%%', uncalibrated_percentages(i)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', ...
         'FontSize', 7, ...
         'Color', 'k');

end


% Add percentage labels above the calibrated bars
%
% NOTE:
% These labels are based on PLACEHOLDER calibrated values and should be
% interpreted only after the actual calibrated growth rates are inserted.

for i = 1:9

    text(x_coords(i,3), ...
         Y(i,3) + (max_y * 0.020), ...
         sprintf('%.1f%%', calibrated_percentages(i)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', ...
         'FontSize', 7, ...
         'FontWeight', 'bold', ...
         'Color', 'k');

end


% Final touches: title and grid

title('Functional Baseline and Calibration Assessment of ecGEMs', ...
      'FontSize', 15, ...
      'FontWeight', 'bold');

grid on;

set(gca, ...
    'GridLineStyle', ':', ...
    'GridAlpha', 0.35, ...
    'Layer', 'top', ...
    'Box', 'on');

hold off;