% setup.m - Add ECOMAP scripts to MATLAB path
%
% Usage:
%   1. Open MATLAB and navigate to the ECOMAP directory
%   2. Run: setup
%   3. Or run from command line: matlab -batch "setup"
%
% This script adds all subdirectories under /scripts to MATLAB's path.
% Run this once per MATLAB session, or add to startup.m for persistence.

fprintf('[setup] Adding ECOMAP scripts to MATLAB path...\n');

% Get the directory where this script is located
scriptDir = fileparts(mfilename('fullpath'));

% Path to scripts directory
scriptsPath = fullfile(scriptDir, 'scripts');

if ~exist(scriptsPath, 'dir')
    error('setup:scripts_not_found', 'Scripts directory not found: %s', scriptsPath);
end

% Add scripts and all subdirectories to path
addpath(genpath(scriptsPath));

% Save path for future sessions (optional)
fprintf('[setup] Saving MATLAB path...\n');
savepath;

fprintf('[setup] Done! Added to path:\n');
fprintf('  %s\n', scriptsPath);
fprintf('\n');
fprintf('NOTE: If you want to permanently add these paths, add the following to your startup.m:\n');
fprintf('  addpath(''%s'', ''-BEGIN'');\n', scriptsPath);
fprintf('  savepath;\n');
