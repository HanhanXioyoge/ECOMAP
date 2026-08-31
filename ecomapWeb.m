function ecomapWeb(varargin)
%ECOMAPWEB Launch the ECOMAP web UI from MATLAB — no console window, no popups.
%
%   ecomapWeb                 % start (if needed) AND open the browser  <- default
%   ecomapWeb('start')        % start AND open the browser
%   ecomapWeb('open')         % open the browser (starts first if needed)
%   ecomapWeb('stop')         % stop the server (kills the process tree)
%   ecomapWeb('status')       % print whether it is running
%   ecomapWeb('restart')      % stop then start
%   ecomapWeb('log')          % print the server log
%
%   Everything is spawned via Java ProcessBuilder with pythonw.exe, so no
%   cmd.exe / console window is ever created on Windows. Your MATLAB session
%   is untouched: the backend runs in its own Python process with its own
%   MATLAB engine.
%
%   Port: 8000 by default, override with setenv('ECOMAP_PORT','8765').
%
%   See also: scripts/web/server/run.py (the backend it starts).

    root = fileparts(mfilename('fullpath'));

    if nargin == 0
        action = 'open';
    elseif islogical(varargin{1})
        if varargin{1}, action = 'open'; else, action = 'start'; end
    else
        action = lower(strtrim(char(varargin{1})));
    end

    switch action
        case 'start',   doStart(root, true);
        case 'open',    doStart(root, true);
        case 'stop',    doStop(root);
        case 'status',  doStatus(root);
        case 'restart', doStop(root); doStart(root, true);
        case 'log',     doLog(root);
        otherwise
            error('ecomapWeb:BadAction', ...
                  'Unknown action "%s". Valid: start, open, stop, status, restart, log.', action);
    end
end

% ---------------------------------------------------------------- actions ---

function doStart(root, alsoOpen)
% See spec §4.1 for the full pipeline.
%
% Pipeline:
%   1. Cross-call protection: refuse to start if a recent failed-marker exists.
%   2. locatePython() — three-priority resolver (config / prompt / auto-search).
%   3. pyExe() — 4-depth venv health check; returns '' if unhealthy.
%   4. bootstrapVenv(root, sysPy) — only when needed; persists hash on success.
%   5. Clear the failed marker on success.
%   6. spawnSilent() the backend.

    % --- 1. Cross-call protection ---
    if recentFailedMarker(root)
        error('ecomapWeb:VenvCorrupt', ...
              ['Previous venv rebuild also failed (within last 5 min). ' ...
               'Investigate and delete %s to retry.\n' ...
               'Common causes:\n' ...
               '  - requirements.txt has bad package versions\n' ...
               '  - pip cannot reach PyPI (network / proxy)\n' ...
               '  - stale .ecomap-python pointing to a Python whose venv cannot be built'], ...
              failedMarkerPath(root));
    end

    p   = webPort();
    url = sprintf('http://127.0.0.1:%d', p);
    if isListening(p)
        fprintf('[ecomapWeb] already serving on %s\n', url);
        if alsoOpen, openBrowser(url); end
        return;
    end

    % --- 2. Resolve system Python ---
    sysPy = locatePython(root);

    % --- 3. Check venv health ---
    py = pyExe(root);
    if isempty(py)
        % --- 4. Rebuild venv ---
        if ~bootstrapVenv(root, sysPy)
            writeFailedMarker(root, sysPy, 'rebuild_failed');
            error('ecomapWeb:NoVenv', 'venv bootstrap failed; see log.');
        end
        py = pyExe(root);
        if isempty(py)
            writeFailedMarker(root, sysPy, 'venv_corrupt');
            error('ecomapWeb:VenvCorrupt', 'venv rebuilt but still unhealthy.');
        end
    end

    % --- 5. Clear failed marker on success ---
    clearFailedMarker(root);

    % --- 6. Spawn backend ---
    server  = fullfile(root, 'scripts', 'web', 'server');
    logF    = fullfile(root, '.web.log');
    pidF    = pidPath(root);
    clearPid(root);
    fprintf('[ecomapWeb] starting backend on %s (MATLAB engine warm-up, ~30s) ...\n', url);
    spawnSilent(py, {fullfile(server, 'run.py')}, server, logF, p, pidF);

    if waitPort(p, 90)
        pid = readPid(root);
        if isempty(pid)
            fprintf('[ecomapWeb] up — %s\n', url);
        else
            fprintf('[ecomapWeb] up — pid=%d, %s\n', pid, url);
        end
        fprintf('[ecomapWeb] log: %s\n', logF);
        if alsoOpen, openBrowser(url); end
    else
        fprintf(2, '[ecomapWeb] backend spawned but port %d never bound.\n', p);
        fprintf(2, '[ecomapWeb] check the log: ecomapWeb(''log'')\n');
    end
end

function doStop(root)
    p   = webPort();
    pid = readPid(root);

    if ~isempty(pid) && isAlive(pid)
        fprintf('[ecomapWeb] stopping pid=%d (and its MATLAB engine) ...\n', pid);
        killTree(pid);
        for k = 1:60
            if ~isAlive(pid), break; end
            pause(0.25);
        end
    elseif ~isempty(pid)
        fprintf('[ecomapWeb] stale pid %d — cleaning up\n', pid);
    else
        fprintf('[ecomapWeb] no managed pid\n');
    end
    clearPid(root);

    if isListening(p)
        fprintf('[ecomapWeb] port %d still busy; freeing it\n', p);
        freePort(p);
        pause(0.5);
        if isListening(p)
            fprintf(2, '[ecomapWeb] WARNING: port %d still busy\n', p);
            return;
        end
    end
    fprintf('[ecomapWeb] stopped\n');
end

function doStatus(root)
    p   = webPort();
    pid = readPid(root);
    if ~isempty(pid) && isAlive(pid)
        if isListening(p), state = 'listening'; else, state = 'starting (no port yet)'; end
        fprintf('[ecomapWeb] running pid=%d, port=%d (%s)\n', pid, p, state);
    elseif isListening(p)
        fprintf('[ecomapWeb] port %d busy but no managed pid (started elsewhere?)\n', p);
    else
        fprintf('[ecomapWeb] not running\n');
    end
end

function doLog(root)
    logF = fullfile(root, '.web.log');
    if exist(logF, 'file') == 2
        fprintf('%s\n', fileread(logF));
    else
        fprintf('[ecomapWeb] no log yet (%s)\n', logF);
    end
end

% ---------------------------------------------------------------- spawning ---

function spawnSilent(exe, argsCell, workDir, logFile, port, pidFile)
%SPAWNSILENT Start a detached child with no console window, log to file.
    requireJvm();
    cmd = javaArray('java.lang.String', 1 + numel(argsCell));
    cmd(1) = java.lang.String(exe);
    for i = 1:numel(argsCell)
        cmd(i + 1) = java.lang.String(argsCell{i});
    end

    pb = java.lang.ProcessBuilder(cmd);
    pb.directory(java.io.File(workDir));
    % Java's env snapshot is taken at JVM start, so ECOMAP_PORT set via setenv()
    % after MATLAB launched would be invisible to the child. Pass it through.
    pb.environment().put('ECOMAP_PORT', num2str(port));
    pb.environment().put('ECOMAP_PID_FILE', pidFile);
    pb.redirectErrorStream(true);
    % MATLAB cannot name the nested class ProcessBuilder$Redirect with dot
    % syntax, so reach its static factories through javaMethod.
    RD = 'java.lang.ProcessBuilder$Redirect';
    pb.redirectOutput(javaMethod('appendTo', RD, java.io.File(logFile)));
    pb.redirectInput(javaMethod('from', RD, java.io.File(nullDevice())));

    % Fire and forget: run.py records its own PID into pidFile, which is how
    % doStop later finds the tree. We deliberately do not hold onto the Java
    % Process handle — the backend must outlive this MATLAB session.
    pb.start();
end

function out = runCapture(exe, argsCell)
%RUNCAPTURE Run a short command, no window, return combined stdout+stderr.
    requireJvm();
    cmd = javaArray('java.lang.String', 1 + numel(argsCell));
    cmd(1) = java.lang.String(exe);
    for i = 1:numel(argsCell)
        cmd(i + 1) = java.lang.String(argsCell{i});
    end
    pb = java.lang.ProcessBuilder(cmd);
    pb.redirectErrorStream(true);
    proc = pb.start();

    rd    = java.io.BufferedReader(java.io.InputStreamReader(proc.getInputStream()));
    parts = {};
    while true
        line = rd.readLine();
        if isempty(line), break; end
        parts{end+1} = char(line); %#ok<AGROW>
    end
    rd.close();
    proc.waitFor();
    out = strjoin(parts, newline);
end

function requireJvm()
    if ~usejava('jvm')
        error('ecomapWeb:NoJvm', ...
              ['ecomapWeb needs the JVM (it spawns the backend windowlessly through it).' ...
               newline 'Start MATLAB without -nojvm, or run the backend directly:' ...
               newline '  scripts/web/server/.venv/Scripts/python scripts/web/server/run.py']);
    end
end

function d = nullDevice()
    if ispc, d = 'NUL'; else, d = '/dev/null'; end
end

function ok = bootstrapVenv(root, sysPy)
%BOOTSTRAPVENV Create or rebuild the Python venv at <root>/scripts/web/server/.venv
%   and install pinned requirements. Best-effort: returns true if the venv is
%   now usable, false (with fprintf) on any failure. Safe to call repeatedly.
%
%   Steps (order is fixed — see spec §4.2 + §0.2):
%     0. kill any stale backend pid tree using this venv (per-PID, never global)
%     1. health check (file + venvRuns + venvHasDeps + requirementsHash) — if
%        everything is healthy, just write the hash (in case it was missing)
%        and return
%     2. remove stale venv (with retry)
%     3. create new venv
%     4. upgrade pip
%     5. install requirements.txt
%     6. venvHasDeps() re-check (must pass)
%     7. write requirements hash
%
%   The caller (doStart) decides whether to write a failed marker.
    ok = false;
    if isempty(sysPy) || exist(sysPy, 'file') ~= 2
        fprintf(2, '[ecomapWeb] bootstrapVenv: sysPy missing or not a file: %s\n', sysPy);
        return;
    end
    venv = fullfile(root, 'scripts', 'web', 'server', '.venv');
    req  = fullfile(root, 'scripts', 'web', 'server', 'requirements.txt');

    % --- Step 0: kill any stale backend pid tree ---
    pid = readPid(root);
    if ~isempty(pid) && isAlive(pid)
        fprintf('[ecomapWeb] killing stale backend pid=%d before rebuilding venv ...\n', pid);
        killTree(pid);
        for k = 1:8
            if ~isAlive(pid), break; end
            pause(0.25);
        end
    end

    % --- Step 1: health check (if already healthy, return early) ---
    pyInVenv = venvPython(venv);
    if exist(pyInVenv, 'file') == 2 && venvHasDeps(pyInVenv) && ~requirementsHashChanged(root)
        fprintf('[ecomapWeb] venv already healthy; nothing to bootstrap.\n');
        writeRequirementsHash(root);
        ok = true;
        return;
    end

    % --- Step 2: remove stale venv ---
    if exist(venv, 'dir') == 7
        fprintf('[ecomapWeb] removing stale venv at %s ...\n', venv);
        if ~removeVenvWithRetry(venv, 3)
            fprintf(2, '[ecomapWeb] cannot remove venv after 3 retries. Delete manually: %s\n', venv);
            return;
        end
    end

    % --- Step 3: create venv ---
    fprintf('[ecomapWeb] creating venv at %s ...\n', venv);
    rc = system(sprintf('"%s" -m venv "%s"', sysPy, venv), '-echo');
    if rc ~= 0
        fprintf(2, '[ecomapWeb] venv creation failed (rc=%d)\n', rc);
        return;
    end

    % --- Step 4: upgrade pip ---
    fprintf('[ecomapWeb] upgrading pip ...\n');
    rc = system(sprintf('"%s" -m pip install -U pip', pyInVenv), '-echo');
    if rc ~= 0
        fprintf(2, '[ecomapWeb] pip self-upgrade failed (rc=%d) - continuing\n', rc);
    end

    % --- Step 5: install requirements ---
    if exist(req, 'file') ~= 2
        fprintf(2, '[ecomapWeb] requirements.txt missing at %s\n', req);
        return;
    end
    fprintf('[ecomapWeb] installing backend requirements (this may take a few minutes) ...\n');
    fprintf('[ecomapWeb]   (tail of output on failure: see the rc=N line below)\n');
    rc = system(sprintf('"%s" -m pip install -r "%s"', pyInVenv, req), '-echo');
    if rc ~= 0
        fprintf(2, '[ecomapWeb] pip install failed (rc=%d)\n', rc);
        return;
    end

    % --- Step 6: re-verify venvHasDeps ---
    if ~venvHasDeps(pyInVenv)
        fprintf(2, '[ecomapWeb] post-install health check failed: fastapi/uvicorn missing.\n');
        return;
    end

    % --- Step 7: write requirements hash ---
    writeRequirementsHash(root);
    fprintf('[ecomapWeb] venv bootstrap complete.\n');
    ok = true;
end

function py = locateSystemPython()
%LOCATESYSTEMPYTHON Find a Python interpreter using PATH search, like `where python`.
%   Runs `where python` (Windows) or `which -a python` (Unix) to enumerate ALL
%   python.exe matches on PATH, then picks the first runnable, in-version one.
%   WindowsApps stubs are logged and skipped. Errors out (no `py=''` fallback)
%   if no runnable Python 3.9-3.12 is found.
    cands = runWherePython();
    if isempty(cands)
        error('ecomapWeb:NoPython', ...
              ['No Python interpreter found on PATH.\n' ...
               '  Tried: `where python` (Windows) / `which -a python` (Unix).\n' ...
               '  Fix:\n' ...
               '    1. Install Python 3.9-3.12 (python.org or Microsoft Store), OR\n' ...
               '    2. Make sure python.exe is on PATH, OR\n' ...
               '    3. Write scripts/web/server/.ecomap-python with your path, e.g.:\n' ...
               '         D:/myPython']);
    end

    for i = 1:numel(cands)
        exe = cands{i};
        if isWindowsAppsPath(exe)
            fprintf('[ecomapWeb] skipping WindowsApps stub: %s\n', exe);
            continue;
        end
        if ~existsFile(exe), continue; end
        if versionGateOK(exe, struct('majorMinor', ''))
            py = exe;
            return;
        end
    end

    error('ecomapWeb:NoPython', ...
          ['No runnable Python 3.9-3.12 found on PATH.\n' ...
           '  Candidates seen: %s\n' ...
           '  Fix: install a supported Python (3.9-3.12) and ensure it is on PATH.'], ...
          strjoin(cands, '; '));
end

function cands = runWherePython()
%RUNWHEREPYTHON Run `where python` (Windows) or `which -a python` (Unix).
%   Returns the trimmed list of matching executable paths. Empty list means
%   nothing was found (rc != 0, or the only line was an `INFO:` placeholder
%   that `where` prints even on success in some locales).
    if ispc
        cmd = 'where python 2>nul';
    else
        cmd = 'which -a python 2>/dev/null';
    end
    [rc, out] = system(cmd);
    if rc ~= 0
        cands = {};
        return;
    end
    rawLines = strsplit(strtrim(out), newline);
    cands = {};
    for i = 1:numel(rawLines)
        line = strtrim(rawLines{i});
        if isempty(line),           continue; end
        if startsWith(line, 'INFO:'), continue; end
        cands{end+1} = line;
    end
end

function tf = existsFile(p)
%EXISTSFILE True iff p names an existing file (handles short/long paths).
    tf = exist(p, 'file') == 2;
end

function [tf, info] = venvRuns(pyExePath)
%VENVRUNS Verify the venv interpreter can boot, and parse sys.executable + sys.version.
%   info = struct('executable','','version','','majorMinor','')
    info = struct('executable', '', 'version', '', 'majorMinor', '');
    tf = false;
    cmd = sprintf('"%s" -c "import sys; print(sys.executable); print(sys.version)" 2>nul', pyExePath);
    [rc, out] = system(cmd);
    if rc ~= 0, return; end
    lines = strsplit(strtrim(out), newline);
    if numel(lines) < 2, return; end
    info.executable = strtrim(lines{1});
    info.version    = strtrim(lines{2});
    parts = sscanf(info.version, '%d.%d');
    if numel(parts) >= 2
        info.majorMinor = sprintf('%d.%d', parts(1), parts(2));
    end
    tf = ~isempty(info.executable) && ~isempty(info.majorMinor);
end

function ok = versionGateOK(exe, info)
%VERSIONGATEOK Return true iff the interpreter at `exe` is between MIN and MAX python.
%   When `info` carries a populated majorMinor field the gate is a pure string
%   compare (no extra subprocess). When majorMinor is empty, we fall back to
%   letting python print sys.version_info[:2] — this happens when the caller
%   did not have a venvRuns() result on hand.
    % --- Python version constraints (per spec §3.2) ---
if nargin < 2
    info.majorMinor = '';
end
v = '';
if isfield(info,'majorMinor')
    v = info.majorMinor;
end
if isempty(v)
    [rc,out] = system(sprintf('"%s" -c "import sys; print(f''{sys.version_info.major}.{sys.version_info.minor}'')"',exe));
    if rc~=0
        ok=false;
        return
    end
    v = strtrim(out);
end
nums = sscanf(v,'%d.%d');
if numel(nums) < 2
    ok = false;
    return;
end
major = nums(1);
minor = nums(2);
ok = (major == 3) && (minor >= 9) && (minor <= 12);
end

function tf = le_to(a, b)
%LE_TO String comparison a <= b for dotted versions like '3.12' vs '3.12'.
    tf = strcmp(a, b) || lt_version(a, b);
end

function tf = ge_to(a, b)
    tf = strcmp(a, b) || lt_version(b, a);
end

function tf = lt_version(a, b)
%LT_VERSION True iff dotted version a < b.
    pa = sscanf(a, '%d.%d');
    pb = sscanf(b, '%d.%d');
    if numel(pa) < 2 || numel(pb) < 2, tf = false; return; end
    if pa(1) ~= pb(1), tf = pa(1) < pb(1); return; end
    tf = pa(2) < pb(2);
end

function exe = locateExe(dirGuess, name)
%LOCATEEXE Convert a Python install directory into the absolute path of
% the interpreter executable (python.exe on Windows, python elsewhere).
    if ispc
        exe = fullfile(dirGuess, 'python.exe');
    else
        exe = fullfile(dirGuess, 'python');
    end
    if ~existsFile(exe)
        % Fallback: trust the launcher name and hope PATH holds.
        exe = name;
    end
end

function tf = isWindowsAppsPath(p)
%ISWINDOWSAPPSPATH True iff p is a Microsoft Store Python stub path.
%   These stubs (under ...\WindowsApps\) launch the Store app instead of
%   running Python. We skip them in the primary search and only fall back
%   to them as a last resort.
    tf = ~isempty(p) && ~isempty(strfind(lower(p), 'windowsapps'));
end

function py = parsePathInput(s)
%PARSEPATHINPUT Normalize user input from interactive prompt to an absolute python.exe path.
%   - trims whitespace and surrounding quotes
%   - normalizes forward slashes to backslashes on Windows
%   - if path ends in 'python' or 'python.exe', treat as full exe path
%   - otherwise treat as a directory and append the platform exe name
%   - returns '' on empty / whitespace-only input
    s = strtrim(s);
    if isempty(s), py = ''; return; end
    % strip surrounding quotes
    if (s(1) == '"' && s(end) == '"') || (s(1) == '''' && s(end) == '''')
        s = s(2:end-1);
        s = strtrim(s);
    end
    if isempty(s), py = ''; return; end
    if ispc
        s = strrep(s, '/', '\');
    end
    [~, name, ext] = fileparts(s);
    if strcmpi(name, 'python') || strcmpi([name ext], 'python.exe')
        py = s;
    else
        if ispc
            py = fullfile(s, 'python.exe');
        else
            py = fullfile(s, 'python');
        end
    end
end

function lines = readTextLines(f)
%READTEXTLINES Read a text file, return trimmed non-empty non-comment lines.
%   Returns {} if the file is missing.
    lines = {};
    if exist(f, 'file') ~= 2, return; end
    raw = fileread(f);
    rawLines = strsplit(strtrim(raw), newline);
    for i = 1:numel(rawLines)
        line = strtrim(rawLines{i});
        if isempty(line) || startsWith(line, '#'), continue; end
        lines{end+1} = line; %#ok<AGROW>
    end
end

function info = emptyPythonConfig()
%EMPTYPYTHONCONFIG Return an empty Python config struct.
    info = struct('python', '', 'version', '', 'type', '');
end

function tf = batchMode()
%BATCHMODE True if MATLAB is running headless (no desktop or MW_BATCH set).
    tf = ~usejava('desktop') || ~isempty(getenv('MW_BATCH'));
end

function p = venvPython(venv)
%VENVPYTHON Absolute path to the python executable inside `venv`.
    if ispc
        p = fullfile(venv, 'Scripts', 'python.exe');
    else
        p = fullfile(venv, 'bin', 'python');
    end
end

function h = sha256OfString(s)
%SHA256OFSTRING Lowercase hex SHA256 of the UTF-8 bytes of s.
    md = java.security.MessageDigest.getInstance('SHA-256');
    bytes = md.digest(uint8(s));
    h = lower(char(reshape(dec2hex(bytes(:), 2)', 1, [])));
end

function normalized = normalizeRequirements(reqFile)
%NORMALIZEREQUIREMENTS Strip blank lines and # comments, trim each line.
%   Used by both requirementsHashChanged() and writeRequirementsHash()
%   so that pure-format changes (extra blank lines, comment edits) do
%   NOT trigger a venv rebuild.
    normalized = '';
    if exist(reqFile, 'file') ~= 2, return; end
    lines = readTextLines(reqFile);
    normalized = strjoin(lines, newline);
end

function tf = venvHasDeps(pyExePath)
%VENVHASDEPS True iff the venv has fastapi and uvicorn installed.
%   NOTE: matlabengine is intentionally NOT checked here. matlabengine is
%   installed differently per MATLAB version and consumed by matlab_bridge.py
%   at server startup. Checking it here would cause false-positive rebuilds
%   in environments where matlab_bridge is configured by a different mechanism.
    cmd = sprintf('"%s" -c "import fastapi, uvicorn" 2>nul', pyExePath);
    [rc, ~] = system(cmd);
    tf = (rc == 0);
end

function f = failedMarkerPath(root)
%FAILEDMARKERPATH Canonical location of the cross-call rebuild-failure marker.
    f = fullfile(root, '.ecomap-rebuild-failed');
end

function tf = recentFailedMarker(root)
%RECENTFAILEDMARKER True iff a failure marker exists and is less than 5 minutes old.
    tf = false;
    f = failedMarkerPath(root);
    if exist(f, 'file') ~= 2, return; end
    lines = readTextLines(f);
    for i = 1:numel(lines)
        if startsWith(lines{i}, 'timestamp=')
            stampStr = strtrim(lines{i}(length('timestamp=')+1:end));
            stamp = datenum(stampStr, 'yyyy-mm-ddTHH:MM:SS');
            if isnan(stamp), return; end
            tf = (now - stamp) < (5 / (24 * 60));
            return;
        end
    end
end

function writeFailedMarker(root, sysPy, reason)
%WRITEFAILEDMARKER Record a cross-call failure for next-call protection.
%   Format (key=value, user-readable, easy to email to support):
%     timestamp=2026-08-01T10:30:00
%     python=C:\xxx\python.exe
%     reason=pip_install_failed
    f = failedMarkerPath(root);
    fid = fopen(f, 'w');
    if fid < 0, return; end
    fprintf(fid, 'timestamp=%s\n', datestr(now, 'yyyy-mm-ddTHH:MM:SS'));
    fprintf(fid, 'python=%s\n', sysPy);
    fprintf(fid, 'reason=%s\n', reason);
    fclose(fid);
end

function clearFailedMarker(root)
%CLEARFAILEDMARKER Delete the failure marker if present; no-op otherwise.
    f = failedMarkerPath(root);
    if exist(f, 'file') == 2
        delete(f);
    end
end

function tf = requirementsHashChanged(root)
%REQUIREMENTSHASHCHANGED True iff current requirements.txt (normalized) hash differs from stored.
    tf = false;
    reqFile  = fullfile(root, 'scripts', 'web', 'server', 'requirements.txt');
    hashFile = fullfile(root, 'scripts', 'web', 'server', '.venv', '.ecomap-requirements-hash');
    if exist(reqFile, 'file') ~= 2 || exist(hashFile, 'file') ~= 2
        return;   % one is missing → cannot compare, do not trigger rebuild
    end
    currentHash = sha256OfString(normalizeRequirements(reqFile));
    storedHash  = strtrim(fileread(hashFile));
    tf = ~strcmpi(currentHash, storedHash);
end

function writeRequirementsHash(root)
%WRITEREQUIREMENTSHASH Persist the SHA256 of the normalized requirements.
    reqFile  = fullfile(root, 'scripts', 'web', 'server', 'requirements.txt');
    hashFile = fullfile(root, 'scripts', 'web', 'server', '.venv', '.ecomap-requirements-hash');
    if exist(reqFile, 'file') ~= 2, return; end
    h = sha256OfString(normalizeRequirements(reqFile));
    fid = fopen(hashFile, 'w');
    if fid < 0, return; end
    fwrite(fid, h);
    fclose(fid);
end

function ok = removeVenvWithRetry(venv, maxTries)
%REMOVEVENWITHRETRY rmdir with retry. Windows file locks (lingering python.exe
%   children, antivirus scans) sometimes cause rmdir to fail; a few retries with
%   a 1s pause between usually succeed.
    ok = false;
    for attempt = 1:maxTries
        if ispc
            [rc, ~] = system(sprintf('rmdir /S /Q "%s" 2>nul', venv));
        else
            [rc, ~] = system(sprintf('rm -rf "%s" 2>/dev/null', venv));
        end
        if rc == 0 && exist(venv, 'dir') ~= 7
            ok = true;
            return;
        end
        if attempt < maxTries
            pause(1.0);
        end
    end
end

function py = promptPythonPath(root)
%PROMPTPYTHONPATH Interactive ask-the-user fallback for locatePython() step 2.
%   Up to 3 rounds. Each round the user may enter:
%     - a path (directory or full python.exe)
%     - 's' to switch to auto-search (returns '')
%     - 'q' to quit (hard-error)
%   Returns the verified python.exe path, or '' to signal the caller to
%   fall through to locateSystemPython() (auto-search).
    py = '';
    for round = 1:3
        if batchMode()
            return;
        end
        if round == 1
            fprintf('[ecomapWeb] No valid Python config found.\n');
        else
            fprintf('[ecomapWeb] Invalid input (%d/3). Try again.\n', round - 1);
        end
        raw = input('[ecomapWeb] Enter Python directory (absolute path), or ''s'' to auto-search, or ''q'' to quit: ', 's');
        if isempty(raw)
            raw = '';
        end
        raw = strtrim(raw);
        if strcmpi(raw, 'q')
            error('ecomapWeb:UserQuit', 'User chose to quit.');
        end
        if strcmpi(raw, 's')
            return;
        end
        cand = parsePathInput(raw);
        if isempty(cand)
            continue;
        end
        if exist(cand, 'file') ~= 2
            fprintf(2, '[ecomapWeb] not found: %s\n', cand);
            continue;
        end
        [tf, info] = venvRuns(cand);
        if ~tf
            fprintf(2, '[ecomapWeb] cannot run: %s\n', cand);
            continue;
        end
        if ~versionGateOK(cand, info)
            fprintf(2, '[ecomapWeb] Python version %s is not in [3.9, 3.12].\n', ...
                    info.majorMinor);
            continue;
        end
        % success — persist for next time
        if ~isempty(info.version)
            writePythonConfig(root, cand, '');
        end
        py = cand;
        return;
    end
    % 3 strikes → signal "go to auto-search"
    fprintf('[ecomapWeb] 3 invalid inputs; falling back to auto-search.\n');
end

function py = locatePython(root)
%LOCATEPYTHON Three-priority Python interpreter resolver. Hard-fails on no result.
%
%   ① Read .ecomap-python (if exists and points to a runnable, in-range
%      Python, use it).
%   ② Interactive prompt (skipped in batch mode). Up to 3 invalid
%      inputs → falls through to ③.
%   ③ Auto-search PATH + well-known Windows paths. If only a WindowsApps
%      stub is found, use it but print a warning.
%
%   On success, the chosen Python is persisted to .ecomap-python for
%   next time (no-op if it was already there).
    info = readPythonConfig(root);
    if ~isempty(info.python)
        if exist(info.python, 'file') == 2
            [tf, vi] = venvRuns(info.python);
            if tf && versionGateOK(info.python, vi)
                if ~isempty(info.version) && ~isempty(vi.version) && ~strcmp(info.version, vi.version)
                    fprintf('[ecomapWeb] warning: .ecomap-python says version %s but found %s; using %s.\n', ...
                            info.version, vi.version, vi.version);
                end
                py = info.python;
                return;
            end
        end
        fprintf('[ecomapWeb] warning: .ecomap-python points to %s but it is missing or unusable; will re-detect.\n', info.python);
    end

    if ~batchMode()
        py = promptPythonPath(root);
        if ~isempty(py)
            return;
        end
    end

    py = locateSystemPython();
    if isempty(py)
        error('ecomapWeb:NoPython', ...
              ['No Python interpreter found.\n' ...
               '  Searched: .ecomap-python, PATH (excluding WindowsApps), well-known paths.\n' ...
               '  Fix:\n' ...
               '    1. Install Python 3.9-3.12 (python.org or Microsoft Store), OR\n' ...
               '    2. Create scripts/web/server/.ecomap-python with your path, e.g.:\n' ...
               '         D:/myPython']);
    end
    writePythonConfig(root, py, '');
    fprintf('[ecomapWeb] auto-detected Python: %s (saved to .ecomap-python)\n', py);
end

function info = readPythonConfig(root)
%READPYTHONCONFIG Read .ecomap-python; supports v1 (bare path) and v2 (key=value).
%   The format is detected by inspecting the first non-comment line: if it
%   contains '=' it is treated as v2 (parsed line by line, unknown keys
%   silently ignored), otherwise the entire file is treated as v1 (the
%   first non-comment line is the python path).
%
%   Returns an emptyPythonConfig() if the file is missing.
    info = emptyPythonConfig();
    configFile = fullfile(root, 'scripts', 'web', 'server', '.ecomap-python');
    lines = readTextLines(configFile);
    if isempty(lines), return; end

    if contains(lines{1}, '=')
        % ---- v2 key=value ----
        for i = 1:numel(lines)
            parts = regexp(lines{i}, '=', 'split', 'once');   % split on first '=' only
            if numel(parts) ~= 2, continue; end
            key   = strtrim(parts{1});
            value = strtrim(parts{2});
            if isfield(info, key)
                info.(key) = value;
            end
            % unknown keys are silently dropped (forward compatibility)
        end
    else
        % ---- v1 bare path ----
        info.python = lines{1};
    end
end

function ok = writePythonConfig(root, pyPath, pyType)
%WRITEPYTHONCONFIG Write a v2-format .ecomap-python file. Returns true on success.
%   The file always contains a `python=<path>` line. If `pyType` is non-empty,
%   a `type=<pyType>` line is added. Existing files are overwritten atomically
%   (write to .tmp, then rename) so a crash mid-write never leaves a half-baked
%   file that would still be valid v1 (worst case is no v2 metadata).
    ok = false;
    configFile = fullfile(root, 'scripts', 'web', 'server', '.ecomap-python');
    if isempty(pyPath), return; end
    if ispc, pyPath = strrep(pyPath, '/', '\'); end

    tmpFile = [configFile '.tmp'];
    fid = fopen(tmpFile, 'w');
    if fid < 0, return; end
    try
        fprintf(fid, '# ECOMAP Python interpreter config (auto-generated)\n');
        fprintf(fid, 'python=%s\n', pyPath);
        if ~isempty(pyType)
            fprintf(fid, 'type=%s\n', pyType);
        end
    catch
        fclose(fid);
        delete(tmpFile);
        return;
    end
    fclose(fid);
    if exist(configFile, 'file') == 2
        delete(configFile);
    end
    ok = (movefile(tmpFile, configFile) == 1);
end

function p = pyExe(root)
%PYEXE Return the venv's windowless interpreter iff the venv is fully healthy.
%   A stale venv (cross-machine copy) → venvRuns() fails → returns ''.
%   A half-broken venv (pip crashed mid-install) → venvHasDeps() fails → ''.
%   A version-mismatched venv (git pull updated requirements.txt) → hash check
%   fails → ''.
%   In all 'no' cases the caller is expected to invoke bootstrapVenv() to
%   rebuild before trying again.
    p = '';
    venv = fullfile(root, 'scripts', 'web', 'server', '.venv');
    if ispc
        cands = {fullfile(venv, 'Scripts', 'pythonw.exe'), ...
                 fullfile(venv, 'Scripts', 'python.exe')};
    else
        cands = {fullfile(venv, 'bin', 'python')};
    end
    for i = 1:numel(cands)
        if exist(cands{i}, 'file') ~= 2, continue; end
        if ~venvRuns(cands{i}),             continue; end   % depth 2
        if ~venvHasDeps(cands{i}),          continue; end   % depth 3
        if requirementsHashChanged(root),   continue; end   % depth 4
        p = cands{i};
        return;
    end
end

% ------------------------------------------------------------ process ctrl ---

function tf = isAlive(pid)
    tf = false;
    if isempty(pid) || pid <= 0, return; end
    if ispc
        out = runCapture('tasklist', {'/FI', sprintf('PID eq %d', pid), '/NH'});
        tf  = contains(out, sprintf('%d', pid));
    else
        out = runCapture('ps', {'-p', sprintf('%d', pid), '-o', 'pid='});
        tf  = ~isempty(strtrim(out));
    end
end

function killTree(pid)
    if ispc
        runCapture('taskkill', {'/F', '/T', '/PID', sprintf('%d', pid)});
    else
        runCapture('kill', {'-TERM', sprintf('-%d', pid)});   % process group
        runCapture('kill', {'-TERM', sprintf('%d', pid)});
    end
end

function freePort(p)
    if ispc
        ps = sprintf(['$c = Get-NetTCPConnection -LocalPort %d -State Listen ' ...
                      '-ErrorAction SilentlyContinue; ' ...
                      'if ($c) { Stop-Process -Id $c.OwningProcess -Force }'], p);
        runCapture('powershell', {'-NoProfile', '-NonInteractive', '-Command', ps});
    else
        runCapture('fuser', {'-k', sprintf('%d/tcp', p)});
    end
end

% ------------------------------------------------------------------- state ---

function p = webPort()
    v = getenv('ECOMAP_PORT');
    if isempty(v)
        p = 8000;
    else
        p = str2double(v);
        if isnan(p) || p <= 0 || p > 65535
            error('ecomapWeb:BadPort', 'ECOMAP_PORT is not a valid port: "%s"', v);
        end
    end
end

function tf = isListening(p)
    tf = false;
    if ~usejava('jvm'), return; end
    sock = [];
    try
        sock = java.net.Socket();
        sock.connect(java.net.InetSocketAddress('127.0.0.1', p), 400);
        tf = true;
    catch
        tf = false;
    end
    if ~isempty(sock)
        try, sock.close(); catch, end %#ok<NOCOM>
    end
end

function tf = waitPort(p, timeoutSec)
    tf    = false;
    tries = ceil(timeoutSec / 0.5);
    for k = 1:tries
        pause(0.5);
        if isListening(p), tf = true; return; end
    end
end

function f = pidPath(root)
%PIDPATH Canonical location of the .web.pid file written by run.py.
    f = fullfile(root, '.web.pid');
end

function pid = readPid(root)
    pid = [];
    f = fullfile(root, '.web.pid');
    if exist(f, 'file') ~= 2, return; end
    v = str2double(strtrim(fileread(f)));
    if ~isnan(v) && v > 0, pid = v; end
end

function clearPid(root)
    f = fullfile(root, '.web.pid');
    if exist(f, 'file') == 2
        delete(f);
    end
end

% ----------------------------------------------------------------- browser ---

function openBrowser(url)
    % 1) Java Desktop — silent, no console, works in MATLAB desktop mode.
    if usejava('jvm')
        try
            dk = java.awt.Desktop.getDesktop();
            dk.browse(java.net.URI(url));
            fprintf('[ecomapWeb] browser opened at %s\n', url);
            return;
        catch
        end
        % 2) rundll32 via ProcessBuilder — GUI subsystem, still no console.
        try
            runCapture('rundll32', {'url.dll,FileProtocolHandler', url});
            fprintf('[ecomapWeb] browser opened at %s\n', url);
            return;
        catch
        end
    end
    % 3) MATLAB's own builtin (this file is ecomapWeb, so `web` is not shadowed).
    try
        web(url, '-browser');
        fprintf('[ecomapWeb] browser opened at %s\n', url);
    catch
        fprintf('[ecomapWeb] open this in your browser: %s\n', url);
    end
end

