function dbStruct = ParseUniProtData(tsvFilePath)

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
    warning('Expected UniProt file was not found.');
end