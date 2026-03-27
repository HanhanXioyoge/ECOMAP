classdef DockerChecker
    % DockerChecker
    % Utility class to check Docker installation and running status
    % for kcat prediction tools (CatPred, DLKcat, UniKP)

    methods (Static)

        function result = checkDocker()
            % checkDocker
            %   Checks if Docker is installed and running
            %
            % Output:
            %   result - struct with fields:
            %       installed   - logical, true if docker command exists
            %       running     - logical, true if docker daemon is running
            %       version     - string, docker version if available
            %       errorMsg    - string, error message if any check failed

            result = struct();
            result.installed = false;
            result.running = false;
            result.version = '';
            result.errorMsg = '';

            % Check 1: docker --version
            try
                [status, cmdout] = system('docker --version');
                if status == 0
                    result.installed = true;
                    % Parse version string
                    vers = regexp(cmdout, 'Docker version ([\d.]+)', 'tokens');
                    if ~isempty(vers)
                        result.version = vers{1}{1};
                    end
                else
                    result.errorMsg = 'Docker command not found. Please install Docker Desktop.';
                    return;
                end
            catch ME
                result.errorMsg = ['Failed to run docker command: ', ME.message];
                return;
            end

            % Check 2: docker ps (tests if daemon is running)
            try
                [status, cmdout] = system('docker ps');
                if status == 0
                    result.running = true;
                else
                    result.errorMsg = 'Docker is installed but not running. Please start Docker Desktop.';
                end
            catch ME
                result.errorMsg = ['Failed to connect to Docker daemon: ', ME.message];
            end
        end

        function ok = waitForDocker(timeout)
            % waitForDocker
            %   Waits for Docker to become ready (up to timeout seconds)
            %
            % Input:
            %   timeout - scalar, maximum wait time in seconds (default: 30)

            if nargin < 1 || isempty(timeout)
                timeout = 30;
            end

            ok = false;
            startTime = tic;

            while toc(startTime) < timeout
                result = DockerChecker.checkDocker();
                if result.installed && result.running
                    ok = true;
                    return;
                end
                pause(2);
            end
        end

        function displayStatus(result, hTextarea)
            % displayStatus
            %   Displays Docker check result in a text area or command window
            %
            % Inputs:
            %   result    - struct from checkDocker()
            %   hTextarea - optional, uitextarea handle to write to

            msg = DockerChecker.getStatusMessage(result);

            if nargin >= 2 && isvalid(hTextarea)
                hTextarea.Value = msg;
            else
                fprintf('%s\n', msg);
            end
        end

        function msg = getStatusMessage(result)
            % getStatusMessage
            %   Formats the check result into a human-readable message

            if result.installed && result.running
                msg = sprintf('Docker Status: Ready\nVersion: %s', result.version);
            elseif result.installed
                msg = sprintf('Docker Warning: Installed but not running\n%s', result.errorMsg);
            else
                msg = sprintf('Docker Error: Not found\n%s', result.errorMsg);
            end
        end

        function openDockerDesktop()
            % openDockerDesktop
            %   Attempts to open Docker Desktop (platform-specific)

            if ismac
                system('open -a Docker');
            elseif isunix
                system('xdg-open docker-desktop://');
            elseif ispc
                system('start docker desktop://');
            else
                error('Cannot detect OS for opening Docker Desktop');
            end
        end
    end
end
