function result=algFseof(model,biomassRxn,targetRxn,iterations,coefficient,outputFile,parameters)
% FSEOF
%   Implements the Flux Scanning based on Enforced Objective Flux algorithm.
%
% Input:
%   model           a model structure
%   biomassRxn      string with reaction ID of the biomass formation or
%                   growth reaction
%   targetRxn       string with reaction ID of target reaction
%   iterations      numeric indicating number of iterations (optional,
%                   default 10)
%   coefficient     numeric indicating ratio of optimal target reaction
%                   flux, must be less than 1 (optional, default 0.9)
%   outputFile      string with output filename (optional, default prints
%                   to command window)
%
% Output:
%   targets         structure with information for identified targets
%       logical     logical array indicating whether a model reaction was
%                   identified as target by FSEOF
%       slope       numeric array with FSEOF slopes for target reactions
%
% This function writes a CSV file or prints to command window.
% If an output has been specified (targets), it will also generate a
% structure indicating for each model reaction whether it is identified by
% FSEOF as a target and the slope of the reaction when switching from
% biomass formation to product formation.
%
% Usage: targets = FSEOF(model, biomassRxn, targetRxn, iterations,...
%                   coefficient, outputFile)
% form raventoolbox
% ------------------- Parameters & inputs -------------------
if nargin < 7 || isempty(parameters)
    parameters = ParameterManager.getParams();
    if isempty(parameters), error('ParameterManager is not set.'); end
end
biomassRxn=char(biomassRxn);
targetRxn=char(targetRxn);

if nargin<4
    iterations=10;
    coefficient=0.9;
end

if nargin <5
    coefficient=0.9;
end

% outputFile default behavior (MDP extension):
%   nargin < 6  (not passed)         -> workspaces/{model}_{targetRxn}/FSEOF/fseof.csv
%   passed ''  (empty char)         -> legacy: print to command window (output = 0)
%   passed a path                   -> write to that path
if nargin < 6
    outputFile = parameters.designDir;
end
if isempty(outputFile)
    output = 0;
else
    output = 1;
end

% --- Re-run cleanup ---
% Replace only the result for this target; other products are preserved.
algoDir = outputFile;
safeTarget = safeTargetName(targetRxn);
csvName = ['fseof_' safeTarget '.csv'];
matName = ['fseof_' safeTarget '_result.mat'];
for oldName = {csvName, matName}
    p = fullfile(algoDir, oldName{1});
    if exist(p, 'file'), delete(p); end
end

%Find out the maximum theoretical yield of target reaction
model=setParam(model,'obj',targetRxn,1);
sol=solveLP(model,1);
targetMax=sol.f*coefficient;   % 90 percent of the theoretical yield

model=setParam(model,'obj',biomassRxn,1);

fseof.results=zeros(length(model.rxns),iterations);
fseof.target=zeros(length(model.rxns),1);
rxnDirection=zeros(length(model.rxns),1);

%Enforce objective flux iteratively
for i=1:iterations
    n=i*targetMax/iterations;
    model=setParam(model,'lb',targetRxn,n);
    
    sol=solveLP(model,1);
    
    fseof.results(:,i)=sol.x;
    
    %Loop through all fluxes and identify the ones that increased upon the
    %enforced objective flux
    for j=1:length(fseof.results)
        if fseof.results(j,1) > 0   %Check the positive fluxes
            
            if i == 1   %The initial round
                rxnDirection(j,1)=1;
                fseof.target(j,1)=1;
            else
                
                if (fseof.results(j,i) > fseof.results(j,i-1)) & fseof.target(j,1)
                    fseof.target(j,1)=1;
                else
                    fseof.target(j,1)=0;
                end
            end
            
        elseif fseof.results(j,1) < 0 %Check the negative fluxes
            
            if i == 1   %The initial round
                rxnDirection(j,1)=-1;
                fseof.target(j,1)=1;
            else
                if (fseof.results(j,i) < fseof.results(j,i-1)) & fseof.target(j,1)
                    fseof.target(j,1)=1;
                else
                    fseof.target(j,1)=0;
                end
            end
            
        end
        
    end
end

%Build result struct. result.rows is filled inside the main loop below
% (one struct per FSEOF-hit reaction, mirroring the printed/written row).
% result.matFile is filled in at save-time below.
result = struct( ...
    'config',     struct('Iterations', iterations, 'Coefficient', coefficient), ...
    'biomassRxn', biomassRxn, ...
    'targetRxn',  targetRxn, ...
    'outputFile', outputFile, ...
    'matFile',    '', ...
    'rows',       struct('slope', {}, 'rowID', {}, 'enzymeID', {}, ...
                         'enzymeName', {}, 'subsystems', {}, ...
                         'direction', {}, 'grRule', {}), ...
    'targets',    struct('logical', [], 'slope', []));

%Generating output (CSV: comma-separated, with header line, header column
% names matching the printed/written order).
formatSpec='%s,%s,%s,%s,%s,%s,%s\n';
headerLine = 'Slope,rowID,Enzyme_ID,Enzyme_Name,Subsystems,Direction,GrRule';
if output == 1    %Output to a file
    outputFile=char(outputFile);
    fid=fopen(fullfile(outputFile, csvName),'w');
    fprintf(fid,'%s\n',headerLine);
else              %Output to screen
    fprintf('%s\n',headerLine);
end

for num=1:length(fseof.target)
    if fseof.target(num,1) == 1
        A0=num2str(abs(fseof.results(num,iterations)-fseof.results(num,1))/abs(targetMax-targetMax/iterations)); %Slope calculation
        A1=num2str(num);                                  %row ID
        A2=char(model.rxns(num));                         %enzyme ID
        A3=char(model.rxnNames(num));                     %enzyme Name
        if isfield(model,'subSystems') && ~isempty(model.subSystems{num});
            if ~any(cellfun(@(x) iscell(x), model.subSystems));
                subSys = cellfun(@(x) {x}, model.subSystems, 'uni', 0);
            else
                subSys = model.subSystems;
            end
            A4=char(strjoin(subSys{num},';'));                   %Subsystems
        else
            A4='';
        end
        A5=num2str(model.rev(num)*rxnDirection(num,1));   %reaction Dirction
        A6=char(model.grRules(num));                      %Gr Rule

        % Capture this row into result.rows (same data as the printed line).
        result.rows(end+1) = struct( ...
            'slope',      str2double(A0), ...
            'rowID',      str2double(A1), ...
            'enzymeID',   A2, ...
            'enzymeName', A3, ...
            'subsystems', A4, ...
            'direction',  str2double(A5), ...
            'grRule',     A6);

        if output == 1    %Output to a file
            fprintf(fid, formatSpec, A0, A1, ...
                    csvEscape(A2), csvEscape(A3), csvEscape(A4), ...
                    A5, csvEscape(A6));
        else              %Output screen
            fprintf(formatSpec, A0, A1, ...
                    csvEscape(A2), csvEscape(A3), csvEscape(A4), ...
                    A5, csvEscape(A6));
        end
    end
end

if output == 1    %Output to a file
    fclose(fid);
end

% Always populate result.targets (was conditional on nargout == 1).
%   .logical mirrors fseof.target as logical
%   .slope   is end-over-end change divided by scan range, for every reaction
result.targets.logical = logical(fseof.target);
result.targets.slope   = abs(fseof.results(:,iterations)-fseof.results(:,1)) ...
                       / abs(targetMax-targetMax/iterations);

% Persist the full result struct to a .mat inside the MDP workspace.
% Always lands in the same workspace dir as the default output file
% (CSV), regardless of where the user pointed `outputFile`.
matDir = outputFile;
if ~isempty(matDir) && ~exist(matDir, 'dir')
    mkdir(matDir);
end
result.matFile = fullfile(matDir, matName);
save(result.matFile, 'result');
fprintf('[fseof] Saved `result` to %s\n', result.matFile);
end


function safe = safeTargetName(targetRxn)
    safe = regexprep(char(targetRxn), '[^A-Za-z0-9_-]+', '_');
    safe = regexprep(safe, '^_+|_+$', '');
    if isempty(safe), safe = 'target'; end
end


function s = csvEscape(s)
% csvEscape  Quote a string for CSV output if it contains a comma, double
%   quote, or newline. Embedded double quotes are doubled per RFC 4180.

    if isempty(s), return; end
    if any(s == ',') || any(s == '"') || any(s == sprintf('\n')) || any(s == sprintf('\r'))
        s = strrep(s, '"', '""');
        s = ['"' s '"'];
    end
end
