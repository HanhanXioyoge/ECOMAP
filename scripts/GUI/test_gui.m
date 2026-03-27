% Test script for ReconstructionApp
addpath(genpath('D:/project/ECOMAP/ECOMAP/scripts/GUI'));

f = figure('Name', 'test', 'Position', [100 100 1100 700], 'Resize', 'off', 'Color', [0.98 0.98 0.98]);
disp('Figure created');

status = uicontrol(f, 'Style', 'text');
disp('Status text created');
status.String = 'Ready';
status.FontSize = 9;
status.ForegroundColor = [0.4 0.4 0.4];
status.Units = 'pixels';
status.Position = [180 20 700 25];
disp('Status text configured');

stepList = uicontrol(f, 'Style', 'listbox');
disp('StepList created');

pause(1);
delete(f);
disp('Test complete');
exit
