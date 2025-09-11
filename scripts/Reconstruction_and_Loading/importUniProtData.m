function dbStruct = importUniProtData(queryInfo, storageDir)
% importUniProtData
%   Downloads and parses UniProt data based on the provided queryInfo
%   structure and stores it at the specified directory.
%
% Inputs:
%   queryInfo     A structure containing the fields:
%                 - uniprot_id
%                 - uniprot_type
%                 - uniprot_geneidfield
%                 - uniprot_reviewed (boolean)
%   storageDir    Path where the downloaded data should be stored.
%
% Output:
%   dbStruct      Structure with field .uniprot containing:
%                 - ID, genes, eccodes, MW, seq

if nargin < 2
    error('Both queryInfo and storageDir must be provided.');
end

% Prepare parameters for UniProt query
query.ID = queryInfo.uniprot_id;
query.category = queryInfo.uniprot_type;
query.geneField = queryInfo.uniprot_geneidfield;

if strcmpi(query.category, 'taxonomy')
    query.category = 'taxonomy_id';
end

if queryInfo.uniprot_reviewed
    reviewFlag = 'reviewed:true+AND+';
else
    reviewFlag = '';
end

% Create storage directory if needed
if ~exist(storageDir, 'dir')
    mkdir(storageDir);
end

% Path to save UniProt TSV data
tsvFilePath = fullfile(storageDir, 'uniprot.tsv');

% Initialize output structure
dbStruct = [];

% Attempt download if data not already present
if ~exist(tsvFilePath, 'file')
    if isempty(query.ID)
        warning('No UniProt ID specified. Skipping download.');
        return;
    end
    fprintf('Downloading UniProt entries for %s %s...\n', query.category, num2str(query.ID));
    downloadURL = ['https://rest.uniprot.org/uniprotkb/stream?query=' ...
                   reviewFlag query.category ':' num2str(query.ID) ...
                   '&fields=accession%2C' query.geneField ...
                   '%2Cec%2Cmass%2Csequence&format=tsv&compressed=false&sort=protein_name%20asc'];
    try
        websave(tsvFilePath, downloadURL, weboptions('Timeout', 30));
        fprintf('Download complete. File saved at: %s\n', tsvFilePath);
    catch
        error(['Failed to download UniProt data. You may try manually using this URL:\n' ...
               downloadURL '\nThen place the file at:\n' tsvFilePath]);
    end
end

% Parse the TSV file
if exist(tsvFilePath, 'file')
    fid = fopen(tsvFilePath, 'r');
    parsedData = textscan(fid, '%q %q %q %q %q', 'Delimiter', '\t', 'HeaderLines', 1);
    fclose(fid);

    dbStruct.ID      = parsedData{1};
    dbStruct.genes   = parsedData{2};
    dbStruct.eccodes = parsedData{3};
    dbStruct.MW      = str2double(parsedData{4});
    dbStruct.seq     = parsedData{5};

    % Check and report duplicates
    [uniqueIDs, uniqueIdx] = unique(dbStruct.ID, 'stable');
    if numel(uniqueIDs) < numel(dbStruct.ID)
        dupIdx = setdiff(1:numel(dbStruct.ID), uniqueIdx);
        fprintf('Found duplicate UniProt entries. Manual inspection may be required.\n');
        disp(dbStruct.ID(dupIdx));
    end
else
    warning('Expected UniProt file was not found after attempted download.');
end
end