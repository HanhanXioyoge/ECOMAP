function UniKP(fileDir, parameters)
% UniKP
% Run kcat prediction inside local Docker image `hanhanxioyoge:kcat_prediction`.
% Input file   : UniKP_input.csv  (in the mounted data folder)
% Output file  : UniKP.csv (written back to the same folder)
%

    % -------------------- Resolve parameters --------------------
    if ~isempty(fileDir)
        dataDir = fileDir;
    elseif nargin < 2 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if ~isfield(parameters,'dataDir') || isempty(parameters.dataDir)
            error('parameters.dataDir is required.');
        else
            dataDir = parameters.dataDir;
            dataDir = fullfile(dataDir, 'kcatData');
        end
        if isempty(parameters)
            error('ParameterManager is not set.');
        end
    end

    dataPath = dataDir;
    outputHost = fullfile(dataPath, 'UniKP.csv');

    % -------------------- macOS PATH fix (optional) -------------
    if ismac
        % include both Intel & Apple Silicon common paths
        setenv('PATH', ['/usr/local/bin' pathsep '/opt/homebrew/bin' pathsep getenv('PATH')]);
    end

    % -------------------- Check Docker availability -------------
    [stDocker, ~] = system('docker --version');
    if stDocker ~= 0
        error(['Docker not found. Please ensure Docker Desktop is installed and running. ', ...
               'On macOS, try launching MATLAB from a terminal so PATH is inherited.']);
    end

    % -------------------- Build and run docker command ----------
    image = 'hanhanxioyoge/kcat_prediction_dlkcat_unikp:v1.0';

    % Container-side fixed POSIX paths
    inputCont  = '/data/UniKP_input.csv';
    outputCont = '/data/UniKP.csv';

    % Compose command: mount host data dir to /data, then call python UniKP.py
    [st,out] = system('nvidia-smi -L');
    gpuOpt = '';
    if st==0 && contains(out,'GPU')
        gpuOpt = '--gpus all ';
    end
    
    cmd = sprintf(['docker run --rm %s -v "%s":/data -w /app/UniKP %s ', ...
                   'python /app/UniKP/UniKP.py %s %s'], ...
                   gpuOpt, dataPath, image, inputCont, outputCont);

    fprintf('Running UniKP in Docker...\n image: %s\n', image);
    [status, cmdout] = system(cmd);

    % -------------------- Check result --------------------------
    if status ~= 0
        error('UniKP failed (exit code %d).\nCommand:\n%s\n\nOutput:\n%s', status, cmd, cmdout);
    end

    if exist(outputHost, 'file') ~= 2
        error('UniKP did not produce output file: %s', outputHost);
    end

    disp('UniKP prediction completed. Output written to UniKP.csv');
end
