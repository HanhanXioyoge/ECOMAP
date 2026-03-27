function launchReconstructionGUI()
    % launchReconstructionGUI
    %   Launch the ECOMAP Reconstruction GUI
    %
    %   This function initializes and opens the graphical user interface
    %   for enzyme-constraint model reconstruction.
    %
    % Usage:
    %   >> launchReconstructionGUI
    %
    % Requirements:
    %   - MATLAB R2022a or later
    %   - COBRA Toolbox (for loadModel, convertecModel)
    %   - Docker Desktop (for kcat predictions)
    %
    % Notes:
    %   - The GUI supports SBML, JSON, and YAML model formats
    %   - For kcat prediction, Docker must be running
    %   - Set ParameterManager before launching if using custom parameters

    % Add GUI folder to path
    guiPath = fileparts(which('launchReconstructionGUI.m'));
    addpath(genpath(guiPath));

    % Check MATLAB version
    v = ver('MATLAB');
    vNum = str2double(regexp(v.Version, '^\d+\.\d+', 'match', 'once'));
    if vNum < 24.1  % R2024a is 24.1
        warning('ECOMAP:Version', ...
            'This GUI is optimized for MATLAB R2024a or later. Some features may not work correctly.');
    end

    % Launch the app
    app = ReconstructionApp();

    % Make figure visible
    app.UIFigure.Visible = 'on';

    fprintf('ECOMAP Reconstruction GUI launched.\n');
    fprintf('Model folder: %s\n', pwd);
end
