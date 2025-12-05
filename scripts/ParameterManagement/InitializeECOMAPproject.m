function InitializeECOMAPproject(name, path)
% InitializeECOMAPproject
%   Function that create a basic project structure for ECOMAP. If a project
%   with the same name exits, the project will not be created
%
% Input:
%   name               an name for the folder struture used in ECOMAP. Also
%                      creates a basic adapter class, which must be manually
%                      adjusted. If not defined, a dialog box will appear.
%   path               a path where to create the folder. If not defined, a
%                      dialog box will appear.
%
% Usage:
%   InitializeECOMAPproject(name, path)

if nargin < 1 || isempty(name)
    prompt = {'Please enter the name of the ECOMAP project (e.g., ecYeastGEM):'};
    dlgtitle = 'Set project name';
    dims = [1 100];
    definput = {'ecModelGEM'};
    opts.Interpreter = 'tex';
    name = inputdlg(prompt,dlgtitle,dims,definput,opts);
    name = char(name);
end

if nargin < 2 || isempty(path)
    path = uigetdir('Project Folder path');
end

fullPath = fullfile(path, name);
% Validate if the project does not exits
if ~exist(fullPath, 'dir')
    % Create the subfolder
    dir = {'data', 'models', 'analysis'};
    for i = 1:length(dir)
        status = mkdir(fullPath, dir{i});
        fid = fopen(fullfile(fullPath, dir{i}, '.keep'), 'w');
        fclose(fid);
    end

    % Read Template.m contents
    fid = fopen(fullfile(findECOMAProot, 'scripts', 'ParameterManagement', 'Template.m'));
    f = fread(fid, '*char')';
    fclose(fid);

    % Replace key values
    f = strrep(f, 'KEY_Template', [name 'ParameterManagement']);
    f = strrep(f, 'KEY_PATH', path);
    f = strrep(f, 'KEY_NAME', name);

    % Save the class file
    filename = fullfile(path, name, [name 'ParameterManagement.m']);
    fid = fopen(filename, 'w');
    fwrite(fid, f);
    fclose(fid);
else
    printOrange('WARNING: A project with the same name exits at the same location. The project was not created.\n');
end
cd(fullPath)
end
