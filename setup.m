function setup()
%SETUP Add all ECOMAP paths to MATLAB.
%
%   Usage:
%     1. Open MATLAB and navigate to the ECOMAP directory.
%     2. Run:  setup
%     3. Or from a shell:  matlab -batch "setup"
%
%   Adds:
%     - scripts and its subfolders
%     - ecYeast, eciML1515, eciCW773, ecHuman (each organism project)
%     - tutorial
%     - tests/matlab (so unit tests can resolve without manual addpath)
%
%   Runtime/cache folders under scripts are filtered out so the MATLAB path
%   does not include Python virtual environments, uploads, or caches.
%
%   Each directory is added only if it actually exists, so this function is
%   safe to call on partial checkouts.

    here = fileparts(mfilename('fullpath'));
    fprintf('[setup] Adding ECOMAP 2.0 paths from %s ...\n', here);

    scriptsDir = fullfile(here, 'scripts');
    addTreeIfPresent(scriptsDir, { ...
        [filesep '.venv'], ...
        [filesep '.uploads'], ...
        [filesep '_uploads'], ...
        [filesep '__pycache__'], ...
        [filesep '.idea'] ...
        });

    % --- four organism projects + tutorials -----------------------------
    addTreeIfPresent(fullfile(here, 'ecYeast'));
    addTreeIfPresent(fullfile(here, 'eciML1515'));
    addTreeIfPresent(fullfile(here, 'eciCW773'));
    addTreeIfPresent(fullfile(here, 'ecHuman'));
    addFolderIfPresent(fullfile(here, 'tutorial'));

    % --- test tree so functiontests() can resolve helpers ---------------
    addTreeIfPresent(fullfile(here, 'tests', 'matlab'));

    fprintf(['[setup] ECOMAP 2.0 paths ready. ' ...
             'Run ecomapWeb(''start'') to launch the UI.\n']);
end

% ---------------------------------------------------------------- helpers ---

function addFolderIfPresent(p)
%ADDFOLDERIFPRESENT add one directory to the MATLAB path if it exists.
    if ~exist(p, 'dir')
        fprintf('[setup]   skip (missing): %s\n', p);
        return;
    end
    addpath(p);
    fprintf('[setup]   added: %s\n', p);
end

function addTreeIfPresent(p, excludedParts)
%ADDTREEIFPRESENT add a directory tree to the MATLAB path if it exists.
    if nargin < 2
        excludedParts = {};
    end
    if ~exist(p, 'dir')
        fprintf('[setup]   skip (missing): %s\n', p);
        return;
    end

    entries = strsplit(genpath(p), pathsep);
    entries = entries(~cellfun('isempty', entries));
    keep = true(size(entries));
    for i = 1:numel(entries)
        keep(i) = ~matchesExcludedPart(entries{i}, excludedParts);
    end
    entries = entries(keep);
    if isempty(entries)
        fprintf('[setup]   skip (empty): %s\n', p);
        return;
    end

    addpath(strjoin(entries, pathsep));
    fprintf('[setup]   added tree: %s\n', p);
end

function tf = matchesExcludedPart(p, excludedParts)
%MATCHESEXCLUDEDPART true when p is inside an ignored runtime/cache folder.
    tf = false;
    for i = 1:numel(excludedParts)
        if contains(p, excludedParts{i})
            tf = true;
            return;
        end
    end
end
