classdef KcatRepo < handle
% KCATREPO
%   A lightweight repository for storing multiple kcat vectors with Umin and Xi values for ONE model.
%   Stored data fields are minimal and stable for later statistics:
%     - modelName : string
%     - organism  : string
%     - groups    : struct array with fields {name, kcat, sluiceParams, note}
%
%   note is defined as: "method | dataTag"
%
%   IMPORTANT: Each kcat set has a corresponding sluice parameters (from GAUKS calibration).
%   The sluiceParams field is a struct with fields:
%     - reactions: cell array of exchange reaction names (e.g., {'EX_glc__D_e', ...})
%     - umin: numeric array of Umin values for each reaction
%     - xi: numeric array of Xi (proteomic cost) values for each reaction
%
%   Recommended usage:
%     repo = KcatRepo("eciML1515_integrated","Escherichia coli");
%     repo.addGroup("Init", kcat0, [], "InitAssign", "DB+DL");
%     repo.addGroup("Bayes_Total", kcat1, sluiceParams, "ABC-SMC", "GrowthRates+UnconstrainedMaxGrowth");
%     repo.save(outDir);
%     repo2 = KcatRepo.loadMat(fullfile(outDir,"eciML1515_integrated_kcatRepo.mat"));

    properties
        modelName (1,1) string = ""
        organism  (1,1) string = ""
        groups                 = struct('name',{},'kcat',{},'sluiceParams',{},'note',{})
    end

    methods
        function obj = KcatRepo(modelName, organism)
            if nargin >= 1, obj.modelName = string(modelName); end
            if nargin >= 2, obj.organism  = string(organism);  end
        end

        function addGroup(obj, name, kcatVec, sluiceParams, method, dataTag, overwrite)
            % ADDGROUP
            %   Add/replace a named kcat vector with corresponding sluice parameters (Umin + Xi).
            %
            % INPUT
            %   name      : group name, e.g. "Init", "Bayes_Total", "GAUKS_afterBayes"
            %   kcatVec   : numeric vector (will be stored as a column)
            %   sluiceParams: struct with fields {reactions, umin, xi}
            %               - reactions: cell array of exchange reaction names
            %               - umin: numeric array of Umin values for each reaction
            %               - xi: numeric array of Xi (proteomic cost) values
            %               Pass empty [] if no sluice params (e.g., initial kcat)
            %   method    : method name for note (string)
            %   dataTag   : data used for note (string)
            %   overwrite : true/false (default false)

            if nargin < 7 || isempty(overwrite), overwrite = false; end
            if nargin < 6, dataTag = ""; end
            if nargin < 5, method  = ""; end
            if nargin < 4, sluiceParams = struct('reactions',{{}}, 'umin',[], 'xi',[]); end

            name   = string(name);
            method = string(method);
            dataTag= string(dataTag);

            if ~isnumeric(kcatVec) || isempty(kcatVec)
                error('KcatRepo:addGroup', 'kcatVec must be a non-empty numeric vector.');
            end
            kcatVec = kcatVec(:); % force column

            % Validate sluiceParams
            if isempty(sluiceParams)
                sluiceParams = struct('reactions',{{}}, 'umin',[], 'xi',[]);
            elseif ~isstruct(sluiceParams) || ~isfield(sluiceParams, 'reactions')
                error('KcatRepo:addGroup', 'sluiceParams must be a struct with fields "reactions", "umin", and "xi".');
            end
            % Ensure xi field exists
            if ~isfield(sluiceParams, 'xi')
                sluiceParams.xi = [];
            end

            % Enforce consistent length across groups (so later stats are easy)
            if ~isempty(obj.groups)
                n0 = numel(obj.groups(1).kcat);
                if numel(kcatVec) ~= n0
                    error('KcatRepo:addGroup', ...
                        'kcat length mismatch: new=%d, existing=%d.', numel(kcatVec), n0);
                end
            end

            note = strtrim(method + " | " + dataTag);

            idx = find(strcmp(string({obj.groups.name}), name), 1);
            if isempty(idx)
                obj.groups(end+1) = struct('name',name,'kcat',kcatVec,'sluiceParams',sluiceParams,'note',note);
            else
                if ~overwrite
                    error('KcatRepo:addGroup', ...
                        'Group "%s" already exists. Use overwrite=true to replace.', name);
                end
                obj.groups(idx).kcat = kcatVec;
                obj.groups(idx).sluiceParams = sluiceParams;
                obj.groups(idx).note = note;
            end
        end

        function [kcatVec, sluiceParams] = getGroup(obj, name)
            % GETGROUP
            %   Retrieve kcat vector and corresponding sluice parameters (Umin + Xi).
            %
            % OUTPUT
            %   kcatVec      : numeric vector
            %   sluiceParams : struct with fields {reactions, umin, xi}
            name = string(name);
            idx = find(strcmp(string({obj.groups.name}), name), 1);
            if isempty(idx)
                error('KcatRepo:getGroup', 'Group "%s" not found.', name);
            end
            kcatVec = obj.groups(idx).kcat;
            sluiceParams = obj.groups(idx).sluiceParams;
        end

        function model = applyGroupToModel(obj, name, baseModel)
            % APPLYGROUPTOMODEL
            %   Apply a saved kcat set and sluice parameters (Umin + Xi) to a model.
            %
            %   INPUT
            %   name      : group name to retrieve
            %   baseModel : model with sluice structure (from applySluiceStructure)
            %
            %   OUTPUT
            %   model     : model with kcat, Umin, and Xi applied
            %
            %   Usage:
            %     repo = KcatRepo.loadMat('eciML1515_integrated_kcatRepo.mat');
            %     model = applySluiceStructure(baseModel, ex_rxn_list);
            %     model = repo.applyGroupToModel('Bayes_FullData', model);

            name = string(name);
            idx = find(strcmp(string({obj.groups.name}), name), 1);
            if isempty(idx)
                error('KcatRepo:applyGroupToModel', 'Group "%s" not found.', name);
            end

            % Apply kcat
            kcatVec = obj.groups(idx).kcat;
            baseModel.enzymeConstraints.kcat = kcatVec;
            model = UpdateSmatrix(baseModel);

            % Apply Umin and Xi if available
            sluiceParams = obj.groups(idx).sluiceParams;
            if ~isempty(sluiceParams) && ~isempty(sluiceParams.reactions)
                for i = 1:length(sluiceParams.reactions)
                    ex_rxn = sluiceParams.reactions{i};
                    umin_val = sluiceParams.umin(i);
                    xi_val = sluiceParams.xi(i);
                    model = setSluiceParams(model, ex_rxn, umin_val, xi_val);
                end
            end
        end

        function [model, appliedName] = applyToModel(obj, baseModel, groupName)
            % APPLYTOMODEL
            %   Apply kcat set(s) and sluice parameters to a model.
            %
            %   INPUT
            %   baseModel  : model with sluice structure (from applySluiceStructure)
            %   groupName  : (optional) specific group name to apply, or 'all'
            %                If not specified or 'all', applies all groups.
            %
            %   OUTPUT
            %   model      : if groupName specified -> model with that group applied
            %               if groupName='all' -> struct with each group applied
            %   appliedName: name of the applied group(s)
            %
            %   Usage:
            %     repo = KcatRepo.loadMat('eciML1515_kcatRepo.mat');
            %
            %     % Apply specific group
            %     model = repo.applyToModel(baseModel, 'Bayes_FullData');
            %
            %     % Apply all groups (returns struct)
            %     models = repo.applyToModel(baseModel, 'all');
            %     model = models.Bayes_FullData;  % access by group name
            %
            %     % List available groups
            %     repo.listGroups();

            if nargin < 3 || isempty(groupName) || strcmpi(groupName, 'all')
                % Apply ALL groups
                if isempty(obj.groups)
                    error('KcatRepo:applyToModel', 'No groups found in repository.');
                end

                groupNames = string({obj.groups.name});
                fprintf('[KcatRepo] Available groups:\n');
                for i = 1:numel(groupNames)
                    fprintf('  [%d] %s\n', i, groupNames(i));
                end
                fprintf('\n');

                % Apply each group and return as struct
                model = struct();
                appliedName = groupNames;
                for i = 1:numel(groupNames)
                    gName = groupNames(i);
                    fprintf('[KcatRepo] Applying group: %s\n', gName);
                    model.(genvarname(gName)) = obj.applyGroupToModel(gName, baseModel);
                end
                fprintf('[KcatRepo] Applied %d groups.\n', numel(groupNames));
            else
                % Apply specific group
                groupName = string(groupName);
                idx = find(strcmp(string({obj.groups.name}), groupName), 1);
                if isempty(idx)
                    error('KcatRepo:applyToModel', ...
                        'Group "%s" not found. Use listGroups() to see available groups.', groupName);
                end
                fprintf('[KcatRepo] Applying group: %s\n', groupName);
                model = obj.applyGroupToModel(groupName, baseModel);
                appliedName = groupName;
                fprintf('[KcatRepo] Done.\n');
            end
        end

        function names = listGroups(obj)
            % LISTGROUPS
            if isempty(obj.groups)
                names = strings(0,1);
            else
                names = string({obj.groups.name}).';
            end
        end

        function s = toStruct(obj)
            % TOSTRUCT
            %   Export ONLY the minimal required fields for robust MAT saving/loading.
            s = struct( ...
                'modelName', obj.modelName, ...
                'organism',  obj.organism, ...
                'groups',    obj.groups ...
            );
        end

        function save(obj, outDir, fileName)
            % SAVE
            %   Save as MAT to outDir. Variable name in file: kcatRepo (struct).
            if nargin < 2 || strlength(string(outDir))==0
                error('KcatRepo:save', 'outDir is required.');
            end
            outDir = char(outDir);
            if ~exist(outDir,'dir'), mkdir(outDir); end

            if nargin < 3 || strlength(string(fileName))==0
                safeModel = regexprep(char(obj.modelName), '[^\w\-]+', '_');
                fileName = [safeModel, '_kcatRepo.mat'];
            end

            kcatRepo = obj.toStruct(); %#ok<NASGU>
            save(fullfile(outDir, fileName), 'kcatRepo', '-v7.3');
        end
    end

    methods(Static)
        function obj = fromStruct(s)
            % FROMSTRUCT
            obj = KcatRepo(s.modelName, s.organism);
            obj.groups = s.groups;
        end

        function obj = loadMat(matFile)
            % LOADMAT
            L = load(matFile, 'kcatRepo');
            if ~isfield(L,'kcatRepo')
                error('KcatRepo:loadMat', 'MAT file does not contain variable "kcatRepo".');
            end
            obj = KcatRepo.fromStruct(L.kcatRepo);
        end
    end
end
