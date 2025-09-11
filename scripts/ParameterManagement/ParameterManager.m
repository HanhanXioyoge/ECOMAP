classdef ParameterManager
    % ParameterManager
    % A simplified static class to load and retrieve parameters from a manager file.
    methods (Static)
        function params = getParams(managerPath)
            %GETPARAMS Load (optional) and return params from manager.
            %   params = getParams()             - returns params from previously set manager
            %   params = getParams(filePath)     - loads manager file, sets it, returns its params
            
            persistent defaultMgr;
            % If a file path is provided, load and cache the manager
            if nargin > 0 && ~isempty(managerPath)
                [folder,name,ext] = fileparts(managerPath);
                if isempty(ext)
                    name = managerPath;
                else
                    if ~strcmp(ext, '.m')
                        error('ParameterManager:getParams', 'Manager must be a .m file.');
                    end
                    if ~isempty(folder)
                        addpath(folder);
                    end
                end
                defaultMgr = feval(name);
            end
            % Ensure a manager is set
            if isempty(defaultMgr)
                error('ParameterManager:getParams', 'No manager file has been loaded.');
            end
            % Return its params field
            if isstruct(defaultMgr) && isfield(defaultMgr, 'params')
                params = defaultMgr.params;
            else
                error('ParameterManager:getParams', 'Loaded manager does not contain a ''params'' field.');
            end
        end
    end
end
