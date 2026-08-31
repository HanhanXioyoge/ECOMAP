function payload = registryPayload(specs)
%REGISTRYPAYLOAD Convert algorithmRegistry entries into a JSON-ready cell array.
%   Shared by mdpAlgorithms and mdpAlgorithmsByTrack. Lives in its own file
%   (not as a local function) so both bridges can reach it -- MATLAB local
%   functions are private to their defining file, which is what left
%   mdpAlgorithmsByTrack un-callable from the engine.
%
%   Returns a 1xN cell of scalar structs. jsonencode turns that into a JSON
%   array. The CALLER wraps it in the standard envelope (see CONTRACT.md);
%   note that make_ok(cellValue) would be wrong -- struct() expands a cell
%   argument into a struct ARRAY -- so callers must assign .result instead.
    payload = cell(1, numel(specs));
    for i = 1:numel(specs)
        s = specs(i);
        params = cell(1, numel(s.params));
        for k = 1:numel(s.params)
            params{k} = s.params(k);
        end
        payload{i} = struct( ...
            'id',         s.id, ...
            'nameKey',    s.nameKey, ...
            'descKey',    s.descKey, ...
            'track',      s.track, ...
            'bridgeName', s.bridgeName, ...
            'supports',   {reshape(s.supports, 1, [])}, ...
            'solverNeed', s.solverNeed, ...
            'runnable',   s.runnable, ...
            'params',     {params});
    end
end
