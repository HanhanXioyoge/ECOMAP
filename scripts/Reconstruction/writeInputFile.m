function writtenTable = writeInputFile(model, kcatList, DeepLearningModel, parameters)
%%
    if nargin < 4 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.');
        end
    end

    % Set the data storage path
    dataDir = parameters.dataDir;

    


