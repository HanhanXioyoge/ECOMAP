function DownloadUniProtData(queryInfo, storageDir)
%   Downloads UniProt data based on the provided queryInfo
%   structure and stores it at the specified directory.
%
% Inputs:
%   queryInfo     A structure containing the fields:
%                 - ID
%                 - type
%                 - geneIDfield
%                 - reviewed (boolean)
%   storageDir    Path where the downloaded data should be stored.
%

% Prepare parameters for UniProt query
query.ID = queryInfo.ID;
query.category = queryInfo.type;
query.geneField = queryInfo.geneIDfield;

if strcmpi(query.category, 'taxonomy')
    query.category = 'taxonomy_id';
end

if queryInfo.reviewed
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
        urlwrite(downloadURL,tsvFilePath,'Timeout',30);
        % websave(tsvFilePath, downloadURL, weboptions('Timeout', 30));
        fprintf('Download complete. File saved at: %s\n', tsvFilePath);
    catch
        error(['Failed to download UniProt data. You may try manually using this URL:\n' ...
               downloadURL '\nThen place the file at:\n' tsvFilePath]);
    end
end
end