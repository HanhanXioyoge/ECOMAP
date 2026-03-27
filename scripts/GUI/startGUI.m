function startGUI()
    % startGUI
    %   Quick launcher for ECOMAP Reconstruction GUI
    %
    % Usage:
    %   startGUI
    %   ecogui          % shorthand alias
    %
    % This is a convenience wrapper around launchReconstructionGUI

    % Get the directory of this function
    funcPath = fileparts(which(mfilename('fullpath')));

    % Add all subfolders to path
    addpath(genpath(funcPath));

    % Launch
    launchReconstructionGUI();
end
