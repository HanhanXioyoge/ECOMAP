function complexInfo = getComplexdata(taxonomicID, parameters)
% getComplexdata  Download or load ComplexPortal data for a given taxonomy
% SYNTAX
%   complexInfo = getComplexdata(taxonomicID, parameters)
% DESCRIPTION
%   If a local cache file 'ComplexPortal.json' does not exist under
%   parameters.reconstructionDir, this function queries the EBI Complex Portal web
%   service (intact/complex-ws) for complexes matching the given
%   taxonomicID, builds a local database and writes it to
%   parameters.reconstructionDir/ComplexPortal.json.
%
%   If the cache file already exists, the function reads it and returns
%   complexInfo as a struct array with the following fields:
%       - complexID      : string (UniProt/Complex AC)
%       - name           : string (complex name)
%       - species        : numeric or struct describing species
%       - geneName       : cell array of gene names (per complex)
%       - protID         : cell array of protein identifiers (per complex)
%       - stochiometry   : numeric array or cell (stoichiometry per subunit)
%       - defined        : numeric flag (0/1/2) as previously stored
%
% INPUTS
%   taxonomicID : numeric or empty (if empty, parameters.taxonomicID is used)
%   parameters  : struct, must contain .reconstructionDir and .taxonomicID (if not passed,
%                 ParameterManager.getParams() is used)
%
% OUTPUT
%   complexInfo : struct array, one element per complex
%
% Note 
% Refer to the GECKO code

    if nargin < 2 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.')
        end
    end

    if nargin<1 || isempty(taxonomicID)
        taxonomicID = parameters.taxonomicID;
    end

    complexdata_Path = fullfile(parameters.reconstructionDir,'ComplexPortal.json');

    if ~isfile(complexdata_Path)
        if isempty(taxonomicID)
            printOrange('WARNING: No taxonomicID specified.');
            return
        elseif taxonomicID == 0
            taxonomicID = [];
        end
        
        webOptions = weboptions('Timeout', 30);
        try
            url1 = 'https://www.ebi.ac.uk/intact/complex-ws/search/*';
            if ~isempty(taxonomicID)
                url1 = [url1 '?facets=species&filters=species:("' num2str(taxonomicID) '")'];
            end
            data = webread(url1,webOptions);
        catch ME
            if (strcmp(ME.identifier,'MATLAB:webservices:HTTP404StatusCodeError'))
                error('Cannot connect to the Complex Portal, perhaps the server is not responding');
            end
        end
        if data.size == 0
            error('No data could be gathered from Complex Portal for the specified taxonomicID.')
        end
        complexData = cell(data.size,7);
        
        for i = 1:data.size
            url2 = 'https://www.ebi.ac.uk/intact/complex-ws/complex/';
            complexID = data.elements(i,1).complexAC;
            try
                temp = webread([url2 complexID],webOptions);
            catch ME
                if (strcmp(ME.identifier,'MATLAB:webservices:HTTP404StatusCodeError'))
                    printOrange(['WARNING: Cannot retrieve the information for ' complexID '.\n']);
                end
                temp = [];
            end
        
            if ~isempty(temp)
                complexData(i,1) = {temp.complexAc};
                complexData(i,2) = {temp.name};
                complexData(i,3) = {temp.species};
        
                idxIntType = find(strcmpi({temp.participants.interactorType}, 'protein'));
        
                % Some complex reported are 'stable complex', then, save the id
                % complex and but set the genes and protein to a empty string.
                if numel(idxIntType) > 0
                    complexData(i,4) = {{temp.participants(idxIntType).name}};
                    complexData(i,5) = {{temp.participants(idxIntType).identifier}};
                else
                    complexData(i,4) = {{temp.participants.name}};
                    complexData(i,5) = {{temp.participants.identifier}};
                end
        
                % Portal complex has two stochiometry values, a minimum and
                % maximum value. Only minimum will be store. In some cases,
                % some complex does not have stochiometry coefficient, then, it
                % will be fill with zeros
                if ~cellfun('isempty',{temp.participants.stochiometry})
                    % For some reason, if there is only one protein in the complex
                    % split function does nor return a cell nx2, instead is 2x1,
                    % then assign an incorrect stochiometry
                    switch numel(idxIntType)
                        case 0 % Contains complexes
                            stochiometry = split({temp.participants.stochiometry}.', ',');
                            complexData(i,7) = {2};
                        case 1 % Contains one protein
                            stochiometry = split({temp.participants(idxIntType).stochiometry}.', ',').';
                            complexData(i,7) = {1};
                        otherwise
                            stochiometry = split({temp.participants(idxIntType).stochiometry}.', ',');
                            complexData(i,7) = {1};
                    end
                    values = str2double(erase(stochiometry(:,1),"minValue: ")).';
                    complexData(i,6) = {values};
                else
                    complexData(i,6) = {repelem(0,numel(complexData{i,4}))};
                    complexData(i,7) = {0};
                end
            end
        end
        fprintf('\n');
        
        % Expand complexes of complexes
        complexComplex = find([complexData{:,7}]==2);
        if ~isempty(complexComplex)
            for i=1:numel(complexComplex)
                subComplex    = complexData{complexComplex(i),5};
                subComplexS   = complexData{complexComplex(i),6};
                subComplexIdx = find(ismember(complexData(:,1),subComplex));
                allGenes = horzcat(complexData{subComplexIdx,4});
                allProts = horzcat(complexData{subComplexIdx,5});
                allStoch = {complexData{subComplexIdx,6}};
                for j=1:numel(subComplex)
                    allStoch{j}=allStoch{j}*subComplexS(j);
                end
                allStoch = horzcat(allStoch{:});
                [allGenes,ia,ic] = unique(allGenes,'stable');
                allProts = allProts(ia);
                allStoch = splitapply(@sum, allStoch', ic);
                complexData{complexComplex(i),4} = allGenes;
                complexData{complexComplex(i),5} = allProts;
                complexData{complexComplex(i),6} = allStoch;
            end
        end
        
        rowHeadings = {'complexID','name','species','geneName','protID','stochiometry','defined'};
        
        complexInfo = cell2struct(complexData, rowHeadings, 2);
        
        % Convert to a JSON file
        jsontxt = jsonencode(cell2table(complexData, 'VariableNames', rowHeadings));
        % Write to a JSON file
        fid = fopen(complexdata_Path, 'w');
        fprintf(fid, '%s', jsontxt);
        fclose(fid);
        fprintf('Model-specific ComplexPortal database stored at %s\n',complexdata_Path);
    else
        % --- Load existing ComplexPortal.json and decode robustly ---
        try
            % Read file as UTF-8 text
            fid = fopen(complexdata_Path, 'r', 'n', 'UTF-8');
            if fid == -1
                error('Could not open %s for reading.', complexdata_Path);
            end
            txt = fread(fid, '*char')';
            fclose(fid);

            if isempty(strtrim(txt))
                error('File %s is empty.', complexdata_Path);
            end

            decoded = jsondecode(txt);
        catch ME
            error('Failed to read or decode JSON file %s:\n%s', complexdata_Path, ME.message);
        end
            complexInfo = decoded;
    end
end
