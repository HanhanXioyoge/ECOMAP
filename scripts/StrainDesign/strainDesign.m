function combined = strainDesign(modelPath, targetRxn, biomassRxn, algoSpec)
% strainDesign  Batch strain design entry point with workspace management.
%   combined = strainDesign(modelPath, targetRxn, biomassRxn, algoSpec)
%   Runs specified algorithms on a model, saves results to workspace.
%
%   Inputs:
%     modelPath  - char, path to model file
%     targetRxn  - char, target reaction ID
%     biomassRxn - char, biomass reaction ID
%     algoSpec   - cell array of struct, each with fields:
%       .algo   - char, algorithm name ('fseof','optknock','optforce',
%                 'oko','okoplus')
%       .params - cell, name-value pairs for the algorithm
%
%   Output:
%     combined - struct with fields: modelName, biomassRxn, targetRxn,
%                timestamp, algorithms
%
%   Workspace structure:
%     workspaces/{safeModelBase}_{safeTargetRxn}/
%       combined_result.mat          <- combined result
%       FSEOF/result.mat             <- per-algorithm result
%       OptKnock/result.mat
%       optForce/result.mat

    arguments
        modelPath  (1,:) char
        targetRxn  (1,:) char
        biomassRxn (1,:) char
        algoSpec   (1,:) cell
    end

    % ---- 1. Load model ----
    model = loadModel(modelPath);
    [~, modelBase] = fileparts(modelPath);

    % ---- 2. Compute workspace directory ----
    wsDir = workspaceDir(modelBase, targetRxn);
    if ~exist(wsDir, 'dir')
        mkdir(wsDir);
    end

    % ---- 3. Load or initialize combined ----
    combinedFile = fullfile(wsDir, 'combined_result.mat');
    if exist(combinedFile, 'file')
        loaded = load(combinedFile, 'combined');
        combined = loaded.combined;
        if ~strcmp(combined.modelName, modelBase)
            warning('strainDesign:ModelMismatch', ...
                    'Workspace was for "%s", now using "%s".', ...
                    combined.modelName, modelBase);
        end
    else
        combined = struct( ...
            'modelName',  modelBase, ...
            'biomassRxn', biomassRxn, ...
            'targetRxn',  targetRxn, ...
            'timestamp',  '', ...
            'algorithms', struct() ...
        );
    end

    % ---- 4. Run each algorithm ----
    for i = 1:numel(algoSpec)
        spec = algoSpec{i};
        algoKey = standardizeAlgoName(spec.algo);
        params  = spec.params;

        switch algoKey
            case 'FSEOF'
                r = algFseof(model, targetRxn, biomassRxn, params{:});
            case 'OptKnock'
                r = algOptknock(model, targetRxn, biomassRxn, params{:});
            case 'optForce'
                r = algOptforce(model, targetRxn, biomassRxn, params{:});
            case 'OKO'
                r = algOko(model, biomassRxn, targetRxn, params{:});
            case 'OKOPlus'
                % First parameter is the predictor-specific interval table/CSV.
                r = algOkoPlus(model, biomassRxn, targetRxn, params{:});
            otherwise
                error('strainDesign:UnknownAlgo', ...
                      'Unknown algorithm: %s', spec.algo);
        end

        % Save per-algorithm result to subfolder
        algoDir = fullfile(wsDir, algoKey);
        if ~exist(algoDir, 'dir')
            mkdir(algoDir);
        end
        algoFile = fullfile(algoDir, 'result.mat');
        save(algoFile, 'r');
        fprintf('[strainDesign] %s saved to %s\n', algoKey, algoFile);

        % Update combined
        combined.algorithms.(algoKey) = r;
        if isfield(r, 'changes')
            fprintf('[strainDesign] %s done: %d changes\n', algoKey, height(r.changes));
        elseif isfield(r, 'candidates') && isfield(r, 'fallbackUsed')
            fprintf('[strainDesign] %s done: %d candidates, fallback=%d\n', ...
                    algoKey, height(r.candidates), r.fallbackUsed);
        else
            fprintf('[strainDesign] %s done\n', algoKey);
        end
    end

    % ---- 5. Save combined ----
    combined.timestamp = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
    save(combinedFile, 'combined');
    fprintf('[strainDesign] Combined saved to %s\n', combinedFile);
end

% === Local functions ===

function d = workspaceDir(modelBase, targetRxn)
% workspaceDir  Compute workspace directory path from model name and target.

    safeModel  = regexprep(modelBase,  '[^\w\-]', '_');
    safeTarget = regexprep(targetRxn,   '[^\w\-]', '_');
    d = fullfile('workspaces', [safeModel '_' safeTarget]);
end

function n = standardizeAlgoName(name)
% standardizeAlgoName  Normalize algorithm name to canonical form.

    switch lower(name)
        case 'fseof',    n = 'FSEOF';
        case 'optknock', n = 'OptKnock';
        case 'optforce', n = 'optForce';
        case 'oko',      n = 'OKO';
        case {'okoplus','oko+'}, n = 'OKOPlus';
        otherwise,       n = name;
    end
end
