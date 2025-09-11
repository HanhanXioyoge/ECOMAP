function [model, noSMILES] = getMetSmiles1(model, parameters)
% FINDMETSMILES Queries PubChem by unique metabolite names to obtain SMILES.
% For any unique metabolite names that remain without a SMILES string,
% the function performs a secondary search using additional annotation data 
% from model.metMiriams. Any newly retrieved SMILES are stored in a local
% SMILES database file (smilesDB.tsv) for future runs.
%
% Input:
%   model - A model structure with fields:
%           - metNames: a cell array of metabolite names.
%           - metMiriams: a cell array of metabolite annotation structures.
%
% Output:
%   model   - The updated model with a new field 'metSmiles' (cell array of SMILES).
%   noSMILES- A cell array of unique metabolite names for which no SMILES could be found.

    % Load parameters once
    if nargin < 2 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.')
        end
    end 

%% STEP 1: Primary Search by Unique Metabolite Names(GECKO's method for obtaining SMILES from pubchem)
    % Clean up metabolite names (remove any unwanted prefixes, e.g. "prot_")
    cleanedNames = regexprep(model.metNames, '^prot_.*', '');
    [uniqueNames, ~, uniqueIdx] = unique(cleanedNames);
    numUnique = numel(uniqueNames);
    uniqueSmiles = repmat({''}, numUnique, 1);

    % Mark empty unique names as matched (no query needed)
    metMatch = false(numUnique, 1);
    metMatch(strcmp(uniqueNames, '')) = true;
    
    % Set up local SMILES database file (in a "data" folder under the current directory)
    smilesDBfile = fullfile(parameters.dataDir,'smilesDB.tsv'); 
    
    % If the local database exists, load its content.
    if exist(smilesDBfile, 'file') == 2
        fID = fopen(smilesDBfile, 'r');
        raw = textscan(fID, '%s %s', 'Delimiter', '\t');
        fclose(fID);
        db.names = raw{1};
        db.smiles = raw{2};
        [metMatch, metIdx] = ismember(uniqueNames, db.names);
        uniqueSmiles(metMatch) = db.smiles(metIdx(metMatch));
        fprintf('Local SMILES database loaded.\n');
    else
        fprintf('Local SMILES database not found.\n');
    end
    
    % Query PubChem for those unique names that have no SMILES yet.
    webOptions = weboptions('Timeout', 30);
    for i = 1:numUnique
        if metMatch(i)
            continue;  % Skip if already found locally
        end
        retry = 0;
        while retry < 10
            try
                % Construct the PubChem REST API URL using the metabolite name.
                url = ['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' uniqueNames{i} '/property/CanonicalSMILES/TXT'];
                result = webread(url, webOptions);
                % Extract the first token if multiple suggestions are returned.
                tokens = regexp(result, '(^\S+)', 'tokens', 'once');
                if ~isempty(tokens)
                    uniqueSmiles{i} = tokens{1};
                else
                    uniqueSmiles{i} = '';
                end
                retry = 15;  % Mark as success
            catch ME
                % If the error is due to a 404, 400, or 500 HTTP status, break out.
                if any(strcmp(ME.identifier, {...
                        'MATLAB:webservices:HTTP404StatusCodeError', ...
                        'MATLAB:webservices:HTTP400StatusCodeError', ...
                        'MATLAB:webservices:HTTP500StatusCodeError'}))
                    uniqueSmiles{i} = '';
                    retry = 15;
                else
                    retry = retry + 1;
                end
            end
            if retry == 10
                error('Cannot reach PubChem. Check your internet connection and try again.');
            end
        end
        % Append the newly retrieved SMILES to the local database.
        fID = fopen(smilesDBfile, 'a');
        fprintf(fID, '%s\t%s\n', uniqueNames{i}, uniqueSmiles{i});
        fclose(fID);
        fprintf('Retrieved SMILES for %s: %s\n', uniqueNames{i}, uniqueSmiles{i});
    end
    
%% STEP 2: Secondary Search Using Additional Annotations
    % For each unique metabolite name that still has an empty SMILES,
    % find its first occurrence in model.metNames, then use that index to access
    % the corresponding annotation in model.metMiriams. For example, if the annotation
    % contains a ChEBI ID, then try to retrieve the SMILES using getSMILESFromChEBI.
    for i = 1:numUnique
        if ~isempty(uniqueSmiles{i})
            continue;  % SMILES already found from primary search.
        end
        % Find all indices in model.metNames that match this unique name.
        idxList = find(strcmp(cleanedNames, uniqueNames{i}));
        if isempty(idxList)
            continue;
        end
        firstIdx = idxList(1);  % Choose the first matching occurrence.
        if firstIdx <= length(model.metMiriams) && ~isempty(model.metMiriams{firstIdx})
            miriam = model.metMiriams{firstIdx};
            % Check if any annotation field contains "chebi" (case-insensitive).
            chebiIdx = find(contains(lower(miriam.name), 'chebi'), 1, 'first');
            if ~isempty(chebiIdx)
                chebiID = miriam.value{chebiIdx};
                chebiID = erase(chebiID, "CHEBI:");
                altSmile = getSMILESFromChEBIUsingInChIKey(chebiID);
                if ~isempty(altSmile)
                    uniqueSmiles{i} = altSmile;
                    fprintf('Alternative search: Retrieved SMILES for %s using ChEBI ID %s.\n', uniqueNames{i}, chebiID);
                end
            end
            % You can add additional alternative search strategies here.
        end
    end
%% Reassemble SMILES for the original metabolite list
    newSmiles = uniqueSmiles(uniqueIdx);
    noSMILES = uniqueNames(cellfun(@isempty, uniqueSmiles));
    successRatio = 1 - (sum(cellfun(@isempty, uniqueSmiles)) / numUnique);
    fprintf('SMILES were found for %.0f%% of unique metabolite names.\n', successRatio*100);

%% Update the model field 'metSmiles'
    if ~isfield(model, 'metSmiles') || all(cellfun(@isempty, model.metSmiles))
        model.metSmiles = newSmiles;
    else
        emptyIndices = cellfun(@isempty, model.metSmiles);
        model.metSmiles(emptyIndices) = newSmiles(emptyIndices);
    end
end

function smiles = getSMILESFromChEBIUsingInChIKey(chebiID)
    % GETSMILESFROMCHEBIUSINGINCHIKEY retrieves the canonical SMILES for a compound
    % given its ChEBI ID. It first retrieves the InChIKey from ChEBI via the libChEBI API,
    % then queries PubChem using that InChIKey.
    %
    % Input:
    %   chebiID - A string representing the ChEBI ID (e.g., 'CHEBI:15377')
    %
    % Output:
    %   smiles  - A string containing the canonical SMILES. Returns an empty string if retrieval fails.
    
    % Step 1: Retrieve the InChIKey from ChEBI using the libChEBI API.
    inchiKey = getInChIKeyFromChEBI(chebiID);
    if isempty(inchiKey)
        warning('No InChIKey found for ChEBI ID %s.', chebiID);
        smiles = '';
        return;
    end
    
    % Step 2: Query PubChem using the retrieved InChIKey to obtain the canonical SMILES.
    smiles = getSMILESFromPubChemByInChIKey(inchiKey);
end

function inchiKey = getInChIKeyFromChEBI(chebiID)
    % GETINCHIKEYFROMCHEBI retrieves the InChIKey for a given ChEBI ID using the libChEBI API.
    %
    % Input:
    %   chebiID - A string representing the ChEBI ID (e.g., 'CHEBI:15377')
    %
    % Output:
    %   inchiKey - A string containing the InChIKey.
    
    % Ensure libChEBI is on the Java classpath using your getChebiEntity function.
    chebiEntity = getChebiEntity(chebiID);
    try
        % Assume the libChEBI ChebiEntity object has a method getInchiKey().
        inchiKeyJava = chebiEntity.getInchiKey();
        inchiKey = char(inchiKeyJava);
    catch ME
        warning('Failed to retrieve InChIKey for ChEBI ID %s: %s', chebiID, ME.message);
        inchiKey = '';
    end
end

function smiles = getSMILESFromPubChemByInChIKey(inchiKey)
    % GETSMILESFROMPUBCHEMBYINCHIKEY retrieves the canonical SMILES for a compound
    % using its InChIKey via the PubChem REST API.
    %
    % Input:
    %   inchiKey - A string representing the InChIKey.
    %
    % Output:
    %   smiles   - A string containing the canonical SMILES, or an empty string if retrieval fails.
    
    % Construct the URL using the InChIKey.
    url = sprintf('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/property/CanonicalSMILES/TXT', inchiKey);
    
    try
        % Set web options to force uncompressed UTF-8 text.
        options = weboptions('ContentType', 'text', 'CharacterEncoding', 'UTF-8', ...
                             'HeaderFields', {'Accept-Encoding', 'identity'});
        response = webread(url, options);
        % If multiple lines are returned, capture the first non-empty token.
        tokens = regexp(response, '(^\S+)', 'tokens', 'once');
        if ~isempty(tokens)
            smiles = tokens{1};
        else
            smiles = '';
        end
    catch ME
        warning('Failed to retrieve SMILES from PubChem for InChIKey %s: %s', inchiKey, ME.message);
        smiles = '';
    end
end
