function model = fillEnzymeInformation(model)
% FILLENZYMESEQUENCES Retrieve enzyme sequences and molecular weights.
%
% This function checks each enzyme in model.enzymeConstraints.enzymes. For any enzyme
% that does not already have a sequence in model.enzymeConstraints.sequence, the function:
%   1. Retrieves the protein sequence and its relative molecular weight using its UniProt ID.
% The sequence is stored in model.enzymeConstraints.sequence and the molecular weight in
% model.enzymeConstraints.mw.
%
% Input:
%   model - A model structure that contains enzyme constraints.
%
% Output:
%   model - The updated model structure with enzyme sequences and molecular weights.

% Ensure the 'sequence' field exists as a cell array.
if ~isfield(model.enzymeConstraints, 'sequence') || isempty(model.enzymeConstraints.sequence)
    model.enzymeConstraints.sequence = repmat({''}, length(model.enzymeConstraints.enzymes), 1);
end

% Ensure the 'mw' field exists as a numeric vector.
if ~isfield(model.enzymeConstraints, 'mw') || isempty(model.enzymeConstraints.mw)
    model.enzymeConstraints.mw = nan(length(model.enzymeConstraints.enzymes), 1);
end

information.id = model.information.uniprotID;
information.type = model.information.uniprot_type;
information.geneIDfield  = model.information.uniprot_type;




% 处理UniProt数据库
uniprotPath = 'uniprot.tsv';
databases.uniprot = [];

% 检查并加载本地数据库
if exist(uniprotPath, 'file')
    fprintf('Loading local UniProt database...\n');
    fid = fopen(uniprotPath, 'r');
    fileContent = textscan(fid, '%q %q %q %q %q', 'Delimiter', '\t', 'HeaderLines', 1);
    fclose(fid);
    
    databases.uniprot.ID = fileContent{1};
    databases.uniprot.seq = fileContent{5};
    databases.uniprot.MW = str2double(fileContent{4});
    
    % 处理重复条目
    [uniqueIDs, uniqueIdx] = unique(databases.uniprot.ID, 'stable');
    if numel(uniqueIDs) < numel(databases.uniprot.ID)
        duplID = setdiff(1:numel(databases.uniprot.ID), uniqueIdx);
        warning('发现重复UniProt条目: %s', strjoin(databases.uniprot.ID(duplID), ', '));
        databases.uniprot.ID = databases.uniprot.ID(uniqueIdx);
        databases.uniprot.seq = databases.uniprot.seq(uniqueIdx);
        databases.uniprot.MW = databases.uniprot.MW(uniqueIdx);
    end
else
    % Try to download the database
    if isfield(model, 'information') && ~isempty(model.information.uniprotID)
        fprintf('Downloading the UniProt database...\n');
        url = ['https://rest.uniprot.org/uniprotkb/stream?query=' uniprotRev ...
                uniprot.type ':' num2str(uniprot.ID) '&fields=accession%2C' uniprot.geneIDfield ...
                '%2Cec%2Cmass%2Csequence&format=tsv&compressed=false&sort=protein_name%20asc'];
        try
            urlwrite(url,uniprotPath,'Timeout',30);
            fprintf('Model-specific UniProt database stored at %s\n',uniprotPath);
        catch
            error(['Download failed, check your internet connection and try again, or manually download: ' url ...
                ' After downloading, store the file as ' uniprotPath])
        end
    end

            
    % 加载下载的数据库
    fid = fopen(uniprotPath, 'r');
    fileContent = textscan(fid, '%q %q %q %q %q', 'Delimiter', '\t', 'HeaderLines', 1);
    fclose(fid);
    
    databases.uniprot.ID = fileContent{1};
    databases.uniprot.seq = fileContent{5};
    databases.uniprot.MW = str2double(fileContent{4});
        catch ME
            warning('数据库下载失败: %s', ME.message);
        end
    end
end

% 填充酶信息
for i = 1:length(model.enzymeConstraints.enzymes)
    enzymeID = model.enzymeConstraints.enzymes{i};
    if isempty(enzymeID) || ~isempty(model.enzymeConstraints.sequence{i})
        continue
    end
    
    % 优先使用本地数据库
    if ~isempty(databases.uniprot)
        idx = find(strcmpi(databases.uniprot.ID, enzymeID), 1);
        if ~isempty(idx)
            model.enzymeConstraints.sequence{i} = databases.uniprot.seq{idx};
            model.enzymeConstraints.mw(i) = databases.uniprot.MW(idx);
            continue
        end
    end
    
    % 在线查询后备方案
    [seq, mw] = getProteinSequence(enzymeID);
    model.enzymeConstraints.sequence{i} = seq;
    model.enzymeConstraints.mw(i) = mw;
end
end

%%
function [seq, mw] = getProteinSequence(uniprotID)
% GETPROTEINSEQUENCE Retrieves the protein sequence and molecular weight for a given UniProt ID.
%
% Input:
%   uniprotID - A string representing the UniProt accession number.
%
% Output:
%   seq - A string containing the protein sequence.
%   mw  - A number representing the relative molecular weight.
%
% This function fetches the data from UniProt using the TSV format to directly obtain both
% the protein sequence and the molecular weight (without calculating the molecular weight from the sequence).

% Construct the URL to fetch the TSV formatted data from UniProt.
url = ['https://rest.uniprot.org/uniprotkb/search?query=accession:', uniprotID, '&format=tsv&fields=sequence,mass'];

try
    % Create web options with ContentType as text, UTF-8 encoding,
    % and force uncompressed response by setting Accept-Encoding to identity.
    options = weboptions('ContentType', 'text', ...
                         'CharacterEncoding', 'UTF-8', ...
                         'HeaderFields', {'Accept-Encoding','identity'});
                     
    % Read the TSV formatted text from the URL.
    tsvText = webread(url, options);
    
    % Split the text into lines.
    lines = strsplit(tsvText, '\n');
    % Remove empty lines.
    lines = lines(~cellfun('isempty', lines));
    
    % Check if there are at least two lines (header and data).
    if length(lines) < 2
        error('No data retrieved from UniProt for %s', uniprotID);
    end
    
    % The second line contains the actual data.
    data = lines{2};
    % Split the data line by tab.
    fields = strsplit(data, '\t');
    if length(fields) < 2
        error('Unexpected data format from UniProt for %s', uniprotID);
    end
    
    % The first field is the protein sequence, the second field is the molecular weight.
    seq = strtrim(fields{1});
    mwStr = strtrim(fields{2});
    mw = str2double(mwStr);
catch ME
    warning('Failed to retrieve protein sequence and MW for %s: %s', uniprotID, ME.message);
    seq = '';
    mw = NaN;
end
end
