function info = setup(projectName)
%SETUP Configure the MATLAB path for ECOMAP.
%
%   SETUP configures the repository root and the MATLAB source folders below
%   scripts/. It deliberately leaves external toolboxes (RAVEN, COBRA,
%   Gurobi, and others) unchanged.
%
%   SETUP(PROJECTNAME) additionally configures one workspace below projects/,
%   for example:
%
%       setup('ecYeast')
%
%   SETUP('all') adds every discovered project. This is convenient for
%   exploration but can make project-local functions with the same name
%   shadow one another, so selecting one project is preferred.
%
%   INFO = SETUP(...) returns the resolved repository root, configured path
%   entries, discovered projects, and the availability of common external
%   solver entry points.
%
%   Only folders containing MATLAB entry points are added. Runtime folders,
%   virtual environments, uploads, caches, model data, and web assets are not
%   placed on the MATLAB path. Re-running SETUP first removes path entries
%   belonging to this checkout, making project switching deterministic.

    if nargin < 1
        projectName = '';
    end
    if ~(ischar(projectName) || (isstring(projectName) && isscalar(projectName)))
        error('ECOMAP:setup:InvalidProject', ...
            'Project name must be a character vector or a string scalar.');
    end
    projectName = strtrim(char(projectName));

    repoRoot = fileparts(mfilename('fullpath'));
    projectsRoot = fullfile(repoRoot, 'projects');
    projectNames = discoverProjects(projectsRoot);

    pathEntries = {repoRoot};
    pathEntries = [pathEntries, collectCodeFolders(fullfile(repoRoot, 'scripts'))];

    testRoot = fullfile(repoRoot, 'tests', 'matlab');
    pathEntries = [pathEntries, collectCodeFolders(testRoot)];

    selectedProjects = resolveProjects(projectName, projectNames);
    for i = 1:numel(selectedProjects)
        projectRoot = fullfile(projectsRoot, selectedProjects{i});
        pathEntries = [pathEntries, collectCodeFolders(projectRoot)]; %#ok<AGROW>
    end

    pathEntries = unique(pathEntries, 'stable');
    removeRepositoryPaths(repoRoot);
    if ~isempty(pathEntries)
        addpath(strjoin(pathEntries, pathsep), '-begin');
    end

    dependencies = dependencyStatus();
    printSummary(repoRoot, pathEntries, projectNames, selectedProjects, dependencies);

    info = struct();
    info.root = repoRoot;
    info.paths = pathEntries;
    info.availableProjects = projectNames;
    info.selectedProjects = selectedProjects;
    info.dependencies = dependencies;
end

function names = discoverProjects(projectsRoot)
%DISCOVERPROJECTS Return project workspaces in stable alphabetical order.
    names = {};
    if ~isfolder(projectsRoot)
        return;
    end

    entries = dir(projectsRoot);
    entries = entries([entries.isdir]);
    for i = 1:numel(entries)
        name = entries(i).name;
        if startsWith(name, '.')
            continue;
        end
        projectRoot = fullfile(projectsRoot, name);
        hasManifest = isfile(fullfile(projectRoot, 'project.json'));
        hasManager = ~isempty(dir(fullfile(projectRoot, '*ParameterManagement.m')));
        if hasManifest || hasManager
            names{end + 1} = name; %#ok<AGROW>
        end
    end

    if ~isempty(names)
        lowerNames = cellfun(@lower, names, 'UniformOutput', false);
        [~, order] = sort(lowerNames);
        names = names(order);
    end
end

function selected = resolveProjects(requested, available)
%RESOLVEPROJECTS Resolve a case-insensitive project name or the value "all".
    if isempty(requested)
        selected = {};
        return;
    end
    if strcmpi(requested, 'all')
        selected = available;
        return;
    end

    match = find(strcmpi(requested, available), 1);
    if isempty(match)
        if isempty(available)
            detail = 'No project workspaces were discovered.';
        else
            detail = sprintf('Available projects: %s.', strjoin(available, ', '));
        end
        error('ECOMAP:setup:UnknownProject', ...
            'Unknown project "%s". %s', requested, detail);
    end
    selected = available(match);
end

function folders = collectCodeFolders(baseFolder)
%COLLECTCODEFOLDERS Find useful MATLAB path entries below BASEFOLDER.
    folders = {};
    if ~isfolder(baseFolder)
        return;
    end

    candidates = strsplit(genpath(baseFolder), pathsep);
    candidates = candidates(~cellfun('isempty', candidates));
    candidates = sortFolders(candidates, baseFolder);

    for i = 1:numel(candidates)
        folder = candidates{i};
        if isExcludedFolder(folder, baseFolder) || ~containsMatlabEntryPoint(folder)
            continue;
        end
        folders{end + 1} = folder; %#ok<AGROW>
    end
end

function folders = sortFolders(folders, baseFolder)
%SORTFOLDERS Put the root first and make descendant ordering deterministic.
    keys = cell(size(folders));
    for i = 1:numel(folders)
        folder = folders{i};
        if pathsEqual(folder, baseFolder)
            keys{i} = '';
        else
            keys{i} = lower(folder);
        end
    end
    [~, order] = sort(keys);
    folders = folders(order);
end

function tf = isExcludedFolder(folder, baseFolder)
%ISEXCLUDEDFOLDER Match complete directory names, not path substrings.
    excludedNames = { ...
        '.git', '.hg', '.svn', '.idea', '.vscode', ...
        '.venv', 'venv', '__pycache__', 'node_modules', ...
        '.uploads', '_uploads', 'uploads', ...
        'cache', 'dist', 'build', 'jobs', 'workspaces'};

    if pathsEqual(folder, baseFolder)
        tf = false;
        return;
    end

    relative = folder(numel(baseFolder) + 1:end);
    relative = regexprep(relative, '^[\\/]+', '');
    parts = regexp(relative, '[\\/]+', 'split');
    lowerParts = cellfun(@lower, parts, 'UniformOutput', false);
    isHidden = cellfun(@(part) startsWith(part, '.'), parts);
    tf = any(isHidden) || any(ismember(lowerParts, excludedNames));
end

function tf = containsMatlabEntryPoint(folder)
%CONTAINSMATLABENTRYPOINT Check source files and package/class containers.
    patterns = {'*.m', '*.mlx', '*.p', ['*.' mexext]};
    for i = 1:numel(patterns)
        if ~isempty(dir(fullfile(folder, patterns{i})))
            tf = true;
            return;
        end
    end

    children = dir(folder);
    children = children([children.isdir]);
    childNames = {children.name};
    tf = any(startsWith(childNames, '+')) || any(startsWith(childNames, '@'));
end

function removeRepositoryPaths(repoRoot)
%REMOVEREPOSITORYPATHS Remove only entries owned by the current checkout.
    entries = strsplit(path, pathsep);
    owned = false(size(entries));
    prefix = [stripTrailingSeparators(repoRoot) filesep];
    for i = 1:numel(entries)
        entry = stripTrailingSeparators(entries{i});
        owned(i) = pathsEqual(entry, repoRoot) || startsWithPath(entry, prefix);
    end
    entries = entries(owned & ~cellfun('isempty', entries));
    if ~isempty(entries)
        rmpath(strjoin(entries, pathsep));
    end
end

function tf = startsWithPath(value, prefix)
%STARTSWITHPATH Apply the platform's path case-sensitivity rules.
    if ispc
        tf = startsWith(value, prefix, 'IgnoreCase', true);
    else
        tf = startsWith(value, prefix);
    end
end

function tf = pathsEqual(left, right)
%PATHSEQUAL Compare paths without insignificant trailing separators.
    left = stripTrailingSeparators(left);
    right = stripTrailingSeparators(right);
    if ispc
        tf = strcmpi(left, right);
    else
        tf = strcmp(left, right);
    end
end

function value = stripTrailingSeparators(value)
%STRIPTRAILINGSEPARATORS Preserve filesystem roots while normalizing paths.
    value = char(value);
    while numel(value) > 1 && any(value(end) == ['/' '\'])
        if ispc && numel(value) == 3 && value(2) == ':'
            break;
        end
        value(end) = [];
    end
end

function status = dependencyStatus()
%DEPENDENCYSTATUS Report entry points used by ECOMAP without changing them.
    status = struct();
    status.solveLP = functionAvailable('solveLP');
    status.cobra = functionAvailable('optimizeCbModel') && ...
        functionAvailable('changeCobraSolver');
    status.gurobi = functionAvailable('gurobi');
    status.parallelComputingToolbox = license('test', 'Distrib_Computing_Toolbox');
end

function tf = functionAvailable(name)
%FUNCTIONAVAILABLE True for MATLAB, P-code, MEX, and built-in functions.
    tf = any(exist(name, 'file') == [2, 3, 6]) || exist(name, 'builtin') == 5;
end

function printSummary(repoRoot, paths, available, selected, dependencies)
%PRINTSUMMARY Keep setup output compact but actionable.
    fprintf('[ECOMAP] Root: %s\n', repoRoot);
    fprintf('[ECOMAP] Configured %d MATLAB path folders.\n', numel(paths));
    if isempty(selected)
        if isempty(available)
            fprintf('[ECOMAP] No project workspaces were discovered.\n');
        else
            fprintf('[ECOMAP] Projects available: %s\n', strjoin(available, ', '));
            fprintf('[ECOMAP] Select one with setup(''<project-name>'').\n');
        end
    else
        fprintf('[ECOMAP] Selected project(s): %s\n', strjoin(selected, ', '));
    end
    fprintf(['[ECOMAP] External capabilities: solveLP=%s, COBRA=%s, ' ...
        'Gurobi=%s, Parallel=%s.\n'], ...
        stateLabel(dependencies.solveLP), stateLabel(dependencies.cobra), ...
        stateLabel(dependencies.gurobi), ...
        stateLabel(dependencies.parallelComputingToolbox));
    fprintf('[ECOMAP] Ready. Run ecomapWeb(''start'') to launch the UI.\n');
end

function label = stateLabel(tf)
%STATELABEL Convert a logical status to concise output text.
    if tf
        label = 'yes';
    else
        label = 'no';
    end
end
