function [kcat_list, rmse, r] = Benchmark(model, fileName, filePath, parameters)
%
% Input:
%   model - A model structure with fields:
%           - metNames: a cell array of metabolite names.
%           - metMiriams: a cell array of metabolite annotation structures.
%   fileName
%
% Output:
%   model   - The updated model with a new field 'metSmiles' (cell array of SMILES).
%   noSMILES- A cell array of unique metabolite names for which no SMILES could be found.
 
    if nargin < 4 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.')
        end
    end

    if nargin < 3 || isempty(filePath)
        filePath = fullfile(parameters.dataDir,'DLKcat.tsv');
    end
    
    if nargin < 2 || isempty(fileName)
        filePath = fullfile(parameters.dataDir,'DLKcat.tsv');
    end

    % Check if Docker is installed
    [checks.status, checks.out] = system('docker --version');
    if checks.status ~= 0
        error(['Docker not found. Please confirm it is installed and running. ' ...
            'If it is, try starting MATLAB from a terminal so the docker command is visible.'])
    end
    
    disp('Running DLKcat prediction, this may take many minutes, especially the first time.')
    status = system(['docker run --rm -v "' fullfile(params.path,'/data') '":/data ghcr.io/sysbiochalmers/dlkcat-gecko:0.1 /bin/bash -c "python DLKcat.py /data/tempDLKcat.tsv /data/tempDLKcatOutput.tsv"']);
    delete(fullfile(params.path,'/data/tempDLKcat.tsv'));
    
    if status == 0 && exist(fullfile(params.path,'data/tempDLKcatOutput.tsv'))
        movefile(fullfile(params.path,'/data/tempDLKcatOutput.tsv'), filePath);
        disp('DKLcat prediction completed.');
    else    
        error('DLKcat encountered an error or it did not create any output file.')
    end

end
