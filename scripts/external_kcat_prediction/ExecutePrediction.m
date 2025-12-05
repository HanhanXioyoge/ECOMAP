function ExecutePrediction(DeepLearningModel, fileDir, parameters)
% ExecutePrediction
%   Dispatch wrapper that triggers one or more deep-learning kcat
%   prediction pipelines by name.
%
%   Given a list of model tags from {'DLKcat','UniKP','CatPred'}, this
%   function normalizes the input into a 1×N cell array of unique tags,
%   validates them, and then calls the corresponding entry-point
%   functions (DLKcat, UniKP, CatPred) in turn.
%
% INPUT
%   DeepLearningModel : string | char | cellstr
%       - A single tag or a list of tags drawn from:
%           'DLKcat', 'UniKP', 'CatPred'
%       - Examples:
%           'DLKcat'
%           ["DLKcat","UniKP"]
%           {'DLKcat','CatPred','CatPred'}   % duplicates are allowed; they
%                                            % will be removed by this function
%
% OUTPUT
%   (none)  — this function has side effects only (it calls other functions).
%
% SIDE EFFECTS / REQUIREMENTS
%   - The functions DLKcat(), UniKP(), and/or CatPred() must exist on the
%     MATLAB path. Each of those should encapsulate its own configuration,
%     I/O, and error handling.
%
% NOTES
%   - This function uses UNIQUE without the 'stable' flag, which may reorder
%     the execution sequence alphabetically (e.g., {'CatPred','DLKcat','UniKP'}).
%     If you need to preserve the original user-specified order, replace:
%         DeepLearningModel = unique(strtrim(DeepLearningModel(:)'));
%     with:
%         DeepLearningModel = unique(strtrim(DeepLearningModel(:)'),'stable');
%
% EXAMPLES
%   % Run DLKcat only
%   ExecutePrediction('DLKcat');
%
%   % Run DLKcat then UniKP (order may be alphabetized by UNIQUE)
%   ExecutePrediction({'DLKcat','UniKP'});
%
%   % Run all three (duplicates removed)
%   ExecutePrediction(["CatPred","DLKcat","CatPred","UniKP"]);

    if nargin < 2 || ~isempty(fileDir)
        dataDir = fileDir;
    else
        if nargin < 3 || isempty(parameters)
            parameters = ParameterManager.getParams();
            dataDir = parameters.dataDir;
            dataDir = fullfile(dataDir, 'kcatData');
            if isempty(parameters)
                error('ParameterManager is not set.');
            end
        end
    end

    % -------------------- Normalize model list --------------------
    % Accept string scalar/array or character vector and convert them
    % into a cell array of character vectors for uniform handling.
    if isstring(DeepLearningModel) || ischar(DeepLearningModel)
        DeepLearningModel = cellstr(DeepLearningModel);
    end

    % Ensure a row cell array, trim surrounding whitespace, and remove
    % duplicates. Note: UNIQUE here can change the order (see NOTES).
    DeepLearningModel = unique(strtrim(DeepLearningModel(:)'));

    % Allowed tags for dispatch. Any other value is considered invalid.
    validTags = {'DLKcat','UniKP','CatPred'};

    % Validate that every requested tag is among the allowed set.
    % If any tag is invalid, throw a descriptive error.
    if any(~ismember(DeepLearningModel, validTags))
        error('DeepLearningModel must be a subset of {DLKcat, UniKP, CatPred}.');
    end

    % -------------------- Dispatch to pipelines -------------------
    % Iterate over the normalized, validated tag list and call the
    % corresponding entry-point function for each prediction pipeline.
    for i = 1:numel(DeepLearningModel)
        tag = DeepLearningModel{i};
        switch tag
            case 'DLKcat'
                % Entry point for the DLKcat prediction pipeline.
                % Must be available on the MATLAB path.
                DLKcat(dataDir);

            case 'UniKP'
                % Entry point for the UniKP prediction pipeline.
                % Must be available on the MATLAB path.
                UniKP(dataDir);

            case 'CatPred'
                % Entry point for the CatPred prediction pipeline.
                % Must be available on the MATLAB path.
                CatPred(dataDir);
        end
    end
end
