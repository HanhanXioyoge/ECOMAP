function jsonStr = mdpAlgorithms(modelType)
%MDPALGORITHMS Registry filtered by model type, as an enveloped JSON string.
%   The envelope .result carries the array of registry entries whose
%   .supports list contains modelType.
%
%   Returns a JSON STRING (not a struct): Python reaches this through
%   matlab_bridge._call_json, which does json.loads(str(response)) and then
%   validates the envelope before reading .result.
    arguments
        modelType (1,:) char
    end
    try
        specs = algorithmRegistry();
        keep  = arrayfun(@(s) ismember(modelType, s.supports), specs);
        entries = registryPayload(specs(keep));
    catch err
        jsonStr = jsonencode(make_err('err_param_invalid', err.message));
        return;
    end
    % Assign .result rather than make_ok(entries): struct('result', aCell)
    % would expand the cell into a 1xN struct ARRAY instead of one envelope.
    envelope = make_ok([]);
    envelope.result = entries;
    jsonStr = jsonencode(envelope);
end
