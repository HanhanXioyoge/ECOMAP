function jsonStr = mdpAlgorithmsByTrack(track)
%MDPALGORITHMSBYTRACK Registry filtered by track, as an enveloped JSON string.
%   track is one of 'recon', 'calib', 'analysis', 'design', or 'all'.
%
%   This used to be a LOCAL function inside mdpAlgorithms.m, which made it
%   invisible to the MATLAB engine -- matlab_bridge.algorithms_by_track would
%   have thrown "Undefined function" on every call. Same failure mode as FIX-4.
    arguments
        track (1,:) char
    end
    try
        specs = algorithmRegistry();
        if strcmpi(track, 'all')
            keep = true(size(specs));
        else
            keep = arrayfun(@(s) strcmp(s.track, track), specs);
        end
        entries = registryPayload(specs(keep));
    catch err
        jsonStr = jsonencode(make_err('err_param_invalid', err.message));
        return;
    end
    % See mdpAlgorithms: assign .result to avoid struct() expanding the cell.
    envelope = make_ok([]);
    envelope.result = entries;
    jsonStr = jsonencode(envelope);
end
