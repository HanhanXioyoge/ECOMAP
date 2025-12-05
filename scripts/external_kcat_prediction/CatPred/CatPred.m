function CatPred(fileDir, parameters)
% CatPred
% Run kcat prediction inside local Docker image `hanhanxioyoge:kcat_prediction_CatPred`.
% Input file   : CatPred_input.csv  (in the mounted data folder)
% Output file  : CatPred.csv (written back to the same folder)
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
    outputHost = fullfile(dataPath, 'CatPred.csv');

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
    image = 'hanhanxioyoge/kcat_prediction_catpred:v1.0';

    % Container-side fixed POSIX paths
    inputCont      = '/data/CatPred_input.csv';
    checkpoint_dir = '/app/data/pretrained/production/kcat/';

    % Compose command: mount host data dir to /data, then call python CatPred.py
    gpuOpt  = '--gpus all';
    gpuFlag = '--use_gpu';

    cmd = sprintf(['docker run --rm %s -v "%s":/data %s ' , ...
                   'python /app/CatPred.py --parameter kcat --input_file %s %s --checkpoint_dir %s'], ...
                   gpuOpt, dataPath, image, inputCont, gpuFlag, checkpoint_dir);

    fprintf('Running CatPred in Docker...\n image: %s\n', image);
    [status, cmdout] = system(cmd);

    % -------------------- Check result --------------------------
    if status ~= 0
        error('CatPred failed (exit code %d).\nCommand:\n%s\n\nOutput:\n%s', status, cmd, cmdout);
    end

    deletefile_1 = fullfile(dataPath, 'CatPred_input_input.csv');
    deletefile_2 = fullfile(dataPath, 'CatPred_input_input.json');
    outputfile   = fullfile(dataPath, 'CatPred_input_input_output.csv');
    targetFile   = fullfile(dataPath, 'CatPred.csv');

    filesToDelete = {deletefile_1, deletefile_2};
    for k = 1:numel(filesToDelete)
        f = filesToDelete{k};
        if isfile(f)
            try
                delete(f);
                fprintf('Deleted: %s\n', f);
            catch ME
                warning('Failed to delete "%s": %s', f, ME.message);
            end
        else
            warning('File not found (skipped): %s', f);
        end
    end
    
    if ~isfile(outputfile)
        error('Required file not found: %s', outputfile);
    end
    
    if isfile(targetFile)
        try
            delete(targetFile);
            fprintf('Removed existing target to avoid conflict: %s\n', targetFile);
        catch ME
            error('Failed to remove existing target "%s": %s', targetFile, ME.message);
        end
    end
    
    [ok, msg, msgID] = movefile(outputfile, targetFile);
    
    if ~ok
        error('Failed to rename "%s" to "%s": [%s] %s', outputfile, targetFile, msgID, msg);
    else
        fprintf('Renamed to: %s\n', targetFile);
    end

    if exist(outputHost, 'file') ~= 2
        warning('CatPred did not produce output file: %s', outputHost);
    end

    disp('CatPred prediction completed. Output written to CatPred.csv');
end
