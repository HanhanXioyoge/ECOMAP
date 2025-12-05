function savePath = saveModel(modelVar, fileName, parameters)
% saveModel  Save a model variable to a .mat file (easy call style).
%
% USAGE
%   savePath = saveModel(ecModel, 'eciML1515.mat', parameters)
%   % Saves ecModel into <parameters.modelDir>\eciML1515.mat
%   % Inside the MAT file, the variable name will be "eciML1515".
%
% INPUTS
%   modelVar   : the model variable to save (struct or any MATLAB variable)
%   fileName   : target file name, may include ".mat" and optional subfolders
%                e.g. 'eciML1515.mat' or 'subdir/eciML1515.mat'
%   parameters : struct; if empty/missing, uses ParameterManager.getParams()
%                must contain parameters.modelDir (base target directory)
%
% OUTPUT
%   savePath   : absolute path to the saved .mat file
%
% BEHAVIOR
%   - If fileName has no extension, '.mat' is appended automatically.
%   - If fileName contains subfolders but is NOT absolute, they are resolved
%     under parameters.modelDir (e.g., modelDir/sub1/sub2/file.mat).
%   - If fileName is an ABSOLUTE path, parameters.modelDir is ignored.
%   - The variable stored inside the MAT is named after the base file name.
%   - Variables >= ~2GB are saved with -v7.3 automatically.

    % -------- Parameters --------
    if nargin < 3 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end
    if ~isfield(parameters, 'modelDir') || isempty(parameters.modelDir)
        error('parameters.modelDir is required.');
    end
    modelDir = parameters.modelDir;

    % -------- Resolve fileName / extension --------
    if nargin < 2 || isempty(fileName)
        % Fallback to the caller's variable name
        nm = inputname(1);
        if isempty(nm), nm = 'model'; end
        fileName = [nm, '.mat'];
    end
    if ~(ischar(fileName) || isstring(fileName))
        error('fileName must be a char or string.');
    end
    fileName = char(fileName);

    [p, base, ext] = fileparts(fileName);
    if isempty(ext), ext = '.mat'; end
    if ~strcmpi(ext, '.mat'), ext = '.mat'; end
    if isempty(base), base = 'model'; end

    % -------- Decide save directory --------
    if isempty(p)
        % No path provided → use modelDir
        saveDir = modelDir;
    else
        if isAbsolutePath(p)
            % Absolute path provided → honor it, ignore modelDir
            saveDir = p;
        else
            % Relative subpath → place under modelDir
            saveDir = fullfile(modelDir, p);
        end
    end

    % -------- Ensure directory exists --------
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    % -------- Variable name inside MAT (from base file name) --------
    varName = matlab.lang.makeValidName(base);
    if isempty(varName), varName = 'model'; end

    % -------- Prepare data for saving --------
    S = struct();
    S.(varName) = modelVar;

    % Choose MAT version based on size
    info   = whos('modelVar');
    useV73 = ~isempty(info) && info.bytes >= 2.0e9;  % ~2GB threshold

    % -------- Save --------
    savePath = fullfile(saveDir, [base, ext]);
    if useV73
        save(savePath, '-struct', 'S', '-v7.3');
    else
        save(savePath, '-struct', 'S');
    end
end

% ===== Helper: detect absolute path (Windows/macOS/Linux) =====
function tf = isAbsolutePath(p)
    % Windows drive letter or UNC
    tf = ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once')) || startsWith(p, '\\');
    % POSIX absolute
    tf = tf || startsWith(p, '/');
end
