classdef KcatRepo < handle
% KCATREPO
%   A lightweight repository for storing multiple kcat vectors for ONE model.
%   Stored data fields are minimal and stable for later statistics:
%     - modelName : string
%     - organism  : string
%     - groups    : struct array with fields {name, kcat, note}
%
%   note is defined as: "method | dataTag"
%
%   Recommended usage:
%     repo = KcatRepo("eciML1515_integrated","Escherichia coli");
%     repo.addGroup("Init", kcat0, "InitAssign", "DB+DL");
%     repo.addGroup("Bayes_Total", kcat1, "ABC-SMC", "GrowthRates+UnconstrainedMaxGrowth");
%     repo.save(outDir);
%     repo2 = KcatRepo.loadMat(fullfile(outDir,"eciML1515_integrated_kcatRepo.mat"));

    properties
        modelName (1,1) string = ""
        organism  (1,1) string = ""
        groups                 = struct('name',{},'kcat',{},'note',{})
    end

    methods
        function obj = KcatRepo(modelName, organism)
            if nargin >= 1, obj.modelName = string(modelName); end
            if nargin >= 2, obj.organism  = string(organism);  end
        end

        function addGroup(obj, name, kcatVec, method, dataTag, overwrite)
            % ADDGROUP
            %   Add/replace a named kcat vector.
            %
            % INPUT
            %   name      : group name, e.g. "Init", "Bayes_Total", "GKC_afterBayes"
            %   kcatVec   : numeric vector (will be stored as a column)
            %   method    : method name for note (string)
            %   dataTag   : data used for note (string)
            %   overwrite : true/false (default false)

            if nargin < 6 || isempty(overwrite), overwrite = false; end
            if nargin < 5, dataTag = ""; end
            if nargin < 4, method  = ""; end

            name   = string(name);
            method = string(method);
            dataTag= string(dataTag);

            if ~isnumeric(kcatVec) || isempty(kcatVec)
                error('KcatRepo:addGroup', 'kcatVec must be a non-empty numeric vector.');
            end
            kcatVec = kcatVec(:); % force column

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
                obj.groups(end+1) = struct('name',name,'kcat',kcatVec,'note',note);
            else
                if ~overwrite
                    error('KcatRepo:addGroup', ...
                        'Group "%s" already exists. Use overwrite=true to replace.', name);
                end
                obj.groups(idx).kcat = kcatVec;
                obj.groups(idx).note = note;
            end
        end

        function kcatVec = getGroup(obj, name)
            % GETGROUP
            name = string(name);
            idx = find(strcmp(string({obj.groups.name}), name), 1);
            if isempty(idx)
                error('KcatRepo:getGroup', 'Group "%s" not found.', name);
            end
            kcatVec = obj.groups(idx).kcat;
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
