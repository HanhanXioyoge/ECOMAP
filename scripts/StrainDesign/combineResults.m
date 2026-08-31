function combined = combineResults(model, biomassRxn, targetRxn, varargin)
% combineResults  Wrap per-algorithm results into a unified struct and save
%   it to the MDP workspace.
%
%   combined = combineResults(model, biomassRxn, targetRxn, ...
%                             'FSEOF',    fseof_result, ...
%                             'OptKnock', optknock_result, ...
%                             'optForce', optforce_result);
%
%   Merge semantics (re-run friendly):
%     - On first run: builds a fresh struct with the supplied fields.
%     - On re-run: loads the existing combined_result.mat (if any), updates
%       only the algorithm fields you pass in this call, and writes back.
%       Algorithm results from prior calls are preserved.
%     - Identity mismatch (model / biomassRxn / targetRxn) on the existing
%       file triggers a warning + takeover. Update the field to the new value.
%
%   Recognised name-value pairs (algorithm names are case-insensitive):
%       'FSEOF'      stored as combined.FSEOF    (default: fseof_result.rows)
%       'OptKnock'   stored as combined.OptKnock (full optknock result)
%       'optForce'   stored as combined.optForce (full optforce result)
%
%   Option (must come BEFORE the FSEOF key in the call):
%       'FSEOFMode', 'rows' (default) -> combined.FSEOF = fseof_result.rows
%       'FSEOFMode', 'full'           -> combined.FSEOF = fseof_result
%
%   Anything else triggers a warning and is ignored.
%
%   Inputs:
%       model       struct from loadModel(...). Its file path is NOT used;
%                   identity is read from model.id / .name / .modelID /
%                   .modelName (first non-empty wins), fallback 'model'.
%       biomassRxn  char, biomass reaction id
%       targetRxn   char, target reaction id
%
%   Output:
%       combined - struct with fields:
%           .model        string identifier from the model struct
%           .biomassRxn
%           .targetRxn
%           .timestamp    char, ISO-like 'yyyy-MM-ddTHH:mm:ss' (refreshed each call)
%           .FSEOF        per the FSEOFMode setting (kept across runs if not re-passed)
%           .OptKnock     full optknock result (kept across runs if not re-passed)
%           .optForce     full optforce result (kept across runs if not re-passed)
%           .combinedFile path to the saved .mat
%
%   Saved to: workspaces/{safeModelId}_{safeTargetRxn}/combined_result.mat
%   (one level above the per-algorithm FSEOF/ subfolder used by fseof.m)

    arguments
        model       (1,1) struct
        biomassRxn  (1,:) char
        targetRxn   (1,:) char
    end

    arguments (Repeating)
        varargin
    end

    % Pull a string identifier from the model struct, same logic as
    % defaultOutputFile() in fseof.m (no shared utility, by design).
    modelBase = localModelBaseName(model);

    % Compute save path up front; we need it both for "load existing" and
    % for the final save.
    safeModel  = regexprep(modelBase, '[^\w\-]', '_');
    safeTarget = regexprep(targetRxn, '[^\w\-]', '_');
    wsDir = fullfile('workspaces', [safeModel '_' safeTarget]);
    if ~exist(wsDir, 'dir')
        mkdir(wsDir);
    end
    combinedFile = fullfile(wsDir, 'combined_result.mat');

    % --- Load-or-init: merge semantics on re-run ---
    if exist(combinedFile, 'file') == 2
        loaded   = load(combinedFile, 'combined');
        combined = loaded.combined;

        % Identity sanity-check; warn and adopt the new values if mismatched.
        if ~strcmp(combined.model, modelBase)
            warning('combineResults:ModelMismatch', ...
                    'Existing combined was for "%s", now using "%s".', ...
                    combined.model, modelBase);
            combined.model = modelBase;
        end
        if ~strcmp(combined.biomassRxn, biomassRxn)
            warning('combineResults:BiomassMismatch', ...
                    'Existing combined biomassRxn="%s", now using "%s".', ...
                    combined.biomassRxn, biomassRxn);
            combined.biomassRxn = biomassRxn;
        end
        if ~strcmp(combined.targetRxn, targetRxn)
            warning('combineResults:TargetMismatch', ...
                    'Existing combined targetRxn="%s", now using "%s".', ...
                    combined.targetRxn, targetRxn);
            combined.targetRxn = targetRxn;
        end
        % Refresh the path in case the workspace dir changed.
        combined.combinedFile = combinedFile;

        % Ensure all 3 algo slots exist (e.g. older files might miss some).
        for slot = {'FSEOF', 'OptKnock', 'optForce'}
            if ~isfield(combined, slot{1})
                combined.(slot{1}) = [];
            end
        end
    else
        % First run: fresh struct.
        combined = struct( ...
            'model',        modelBase, ...
            'biomassRxn',   biomassRxn, ...
            'targetRxn',    targetRxn, ...
            'timestamp',    '', ...
            'FSEOF',        [], ...
            'OptKnock',     [], ...
            'optForce',     [], ...
            'combinedFile', combinedFile);
    end

    % --- Update only the algorithm fields the user passed in ---
    fseofMode = 'rows';

    i = 1;
    while i <= numel(varargin)
        if i+1 > numel(varargin)
            error('combineResults:OddArgs', ...
                  'varargin pairs are incomplete; expected name, value, ...');
        end
        key = varargin{i};
        val = varargin{i+1};
        switch lower(key)
            case 'fseofmode'
                fseofMode = lower(char(val));
            case 'fseof'
                switch lower(fseofMode)
                    case 'rows'
                        if isstruct(val) && isfield(val, 'rows')
                            combined.FSEOF = val.rows;
                        else
                            combined.FSEOF = val;   % fallback: keep what we got
                        end
                    otherwise   % 'full' or anything else
                        combined.FSEOF = val;
                end
            case 'optknock'
                combined.OptKnock = val;
            case 'optforce'
                combined.optForce = val;
            otherwise
                warning('combineResults:UnknownKey', ...
                        'Ignoring unknown name: %s', key);
        end
        i = i + 2;
    end

    % Refresh timestamp on every call.
    combined.timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));

    save(combinedFile, 'combined');
    fprintf('[combineResults] Saved to %s (FSEOF=%d, OptKnock=%d, optForce=%d)\n', ...
            combinedFile, ...
            ~isempty(combined.FSEOF), ...
            ~isempty(combined.OptKnock), ...
            ~isempty(combined.optForce));
end

% === Local functions ===

function out = localModelBaseName(model)
% localModelBaseName  Pick a string identifier from a model struct.
%   Mirrors defaultOutputFile() in fseof.m (no shared utility, by spec).
%   Looks at id / name / modelID / modelName in order; first non-empty
%   char/string scalar wins. Falls back to 'model'.

    candidate = '';
    for fld = {'id', 'name', 'modelID', 'modelName'}
        if isfield(model, fld{1}) && ~isempty(model.(fld{1}))
            v = model.(fld{1});
            if ischar(v) || (isstring(v) && isscalar(v))
                candidate = char(v);
                break;
            end
        end
    end
    if isempty(candidate)
        candidate = 'model';
    end
    out = candidate;
end
