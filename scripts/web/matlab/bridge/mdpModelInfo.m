function jsonStr = mdpModelInfo(path)
%MDPMODELINFO Load a model and return an enveloped JSON string.
%   The envelope .result carries: modelType + rxn/met/gene counts + rxn ids.
%   Python reads envelope["result"] via matlab_bridge.model_info, so a bare
%   jsonencode(info) here 500s the request (BridgeContractError).
    arguments
        path (1,:) char
    end
    try
        model = loadModel(path);
    catch err
        jsonStr = jsonencode(make_err('err_model_format', err.message));
        return;
    end
    nMets = 0;  if isfield(model,'mets'),  nMets = numel(model.mets);  end
    nGenes = 0; if isfield(model,'genes'), nGenes = numel(model.genes); end
    info = struct( ...
        'modelType', model.modelType, ...
        'nRxns', numel(model.rxns), ...
        'nMets', nMets, ...
        'nGenes', nGenes, ...
        'rxns', {reshape(cellstr(model.rxns),1,[])} );
    jsonStr = jsonencode(make_ok(info));
end
