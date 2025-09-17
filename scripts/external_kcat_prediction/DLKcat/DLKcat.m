function DLKcat(dataPath, parameters)

    if nargin < 2 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.')
        end
    end
    
    if nargin < 1 || isempty(dataPath)
        dataPath = fullfile(parameters.dataDir,'DLKcat_input.csv');
    elseif strcmp(dataPath(end),{'\','/'})
        dataPath = fullfile(dataPath,'DLKcat_input.csv');
    end

    copyfile(dataPath, fullfile(parameters.dataDir,'tempDLKcat_input.csv'));
        
    %% Check and install requirements
    % On macOS, Docker might not be properly loaded if MATLAB is started via
    % launcher and not terminal.
    if ismac
        setenv('PATH', strcat('/usr/local/bin', ':', getenv("PATH")));
    end
    
    % Check if Docker is installed
    [checks.docker.status, checks.docker.out] = system('docker --version');
    if checks.docker.status ~= 0
        error(['Docker not found. Please confirm it is installed and running. ' ...
            'If it is, try starting MATLAB from a terminal so the docker command is visible.'])
    end
    
    disp('Running DLKcat prediction, this may take many minutes, especially the first time.')
    status = system(['docker run --rm -v "' parameters.dataDir '":/data ghcr.io/sysbiochalmers/dlkcat-gecko:0.1 /bin/bash -c "python DLKcat.py /data/tempDLKcat.tsv /data/tempDLKcatOutput.tsv"']);
    delete(fullfile(parameters.dataDir,'tempDLKcat_input.csv'));
    
    if status == 0 && exist(fullfile(parameters.dataDir,'tempDLKcat_input.csv'))
        movefile(fullfile(parameters.dataDir,'tempDLKcat_input.csv'), dataPath);
        disp('DLKcat prediction completed.');
    else    
        error('DLKcat encountered an error or it did not create any output file.')
    end
end
