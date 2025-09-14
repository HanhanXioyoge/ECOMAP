function writtenTable = writeInputFile(model, DeepLearningModel, parameters)
%%
    if nargin < 3 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.');
        end
    end
    
    % Set the data storage path
    dataDir = parameters.dataDir;
    
    if nargin < 2 || isempty(DeepLearningModel)
        error('writeInputFile:MissingDLModelList', ...
              'DeepLearningModel must be a non-empty list. Allowed: {DLKCAT, UNIKP, CATPRED}.');
    end
    
    if ischar(DeepLearningModel) || isstring(DeepLearningModel)
        DeepLearningModel = cellstr(DeepLearningModel);
    elseif ~iscell(DeepLearningModel)
        error('DeepLearningModel must be char/string/string array/cellstr.');
    end
    
    if ~isfield(model, 'enzymeConstraints') || isempty(model.enzymeConstraints)
        error('The model does not contain the enzymeConstraints structure.');
    end
    
    
    {dataDir, DeepLearningModel};
end