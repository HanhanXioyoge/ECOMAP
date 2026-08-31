function InitializeECOMAPproject(name, path)
% InitializeECOMAPproject
%   Create the ECOMAP project workspace under <root>/projects/<name>.
%
% Input:
%   name     Project folder name. The parameter manager function name is
%            sanitized with matlab.lang.makeValidName when needed.
%   path     ECOMAP root directory. If omitted, a dialog box will appear.
%
% Usage:
%   project_name = 'eciCW773';
%   project_path = findECOMAProot;
%   InitializeECOMAPproject(project_name, project_path)

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

if isequal(path, 0)
    error('InitializeECOMAPproject:NoProjectPath', 'No project path selected.');
end

if iscell(name)
    if isempty(name)
        error('InitializeECOMAPproject:InvalidName', 'Project name cannot be empty.');
    end
    name = name{1};
end
name = strtrim(char(name));
if isempty(name)
    error('InitializeECOMAPproject:InvalidName', 'Project name cannot be empty.');
end

rootPath = char(path);
projectsRoot = fullfile(rootPath, 'projects');
fullPath = fullfile(projectsRoot, name);
safeName = matlab.lang.makeValidName(name);
managerFunction = [safeName 'ParameterManagement'];
managerFile = [managerFunction '.m'];

ensureDir(projectsRoot);
ensureDir(fullPath);

topLevelDirs = {'models', 'reconstruction', 'calibration', 'analysis', 'design'};
for i = 1:numel(topLevelDirs)
    ensureDir(fullfile(fullPath, topLevelDirs{i}));
end

templatePath = fullfile(findECOMAProot, 'scripts', 'ParameterManagement', 'Template.m');
if ~exist(templatePath, 'file')
    error('InitializeECOMAPproject:MissingTemplate', 'Template file not found: %s', templatePath);
end

paramFile = fullfile(fullPath, managerFile);
if ~exist(paramFile, 'file')
    f = readText(templatePath);
    f = strrep(f, 'KEY_Template', managerFunction);
    f = strrep(f, 'KEY_PATH', rootPath);
    f = strrep(f, 'KEY_NAME', name);
    writeText(paramFile, f);
end

projectJson = fullfile(fullPath, 'project.json');
if ~exist(projectJson, 'file')
    project = struct();
    project.projectName = name;
    project.projectDir = fullPath;
    project.parameterManager = managerFile;
    project.parameterManagerFunction = managerFunction;
    project.stage = 'Initialized';
    project.createdAt = timestampUTC();
    project.updatedAt = project.createdAt;
    project.directories = struct( ...
        'models', fullfile(fullPath, 'models'), ...
        'reconstruction', fullfile(fullPath, 'reconstruction'), ...
        'calibration', fullfile(fullPath, 'calibration'), ...
        'analysis', fullfile(fullPath, 'analysis'), ...
        'design', fullfile(fullPath, 'design'));
    writeText(projectJson, jsonencode(project));
end

cd(fullPath)
end

function ensureDir(pathToCreate)
if ~exist(pathToCreate, 'dir')
    [ok, msg] = mkdir(pathToCreate);
    if ~ok
        error('InitializeECOMAPproject:MkdirFailed', 'Could not create %s: %s', pathToCreate, msg);
    end
end
end

function txt = readText(filename)
fid = fopen(filename, 'r');
if fid < 0
    error('InitializeECOMAPproject:ReadFailed', 'Could not read %s.', filename);
end
cleanup = onCleanup(@() fclose(fid));
txt = fread(fid, '*char')';
delete(cleanup);
end

function writeText(filename, txt)
fid = fopen(filename, 'w');
if fid < 0
    error('InitializeECOMAPproject:WriteFailed', 'Could not write %s.', filename);
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, txt);
delete(cleanup);
end

function stamp = timestampUTC()
try
    stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
catch
    stamp = datestr(now, 'yyyy-mm-ddTHH:MM:SSZ');
end
end
