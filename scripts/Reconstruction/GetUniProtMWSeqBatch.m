function [foundIds, MWs, seqs, foundIdx] = GetUniProtMWSeqBatch(accessions, varargin)
% GetUniProtMWSeqBatch  Fetch only entries that have BOTH mass and sequence.
%
% [foundIds, MWs, seqs, foundIdx] = GetUniProtMWSeqBatch(accessions)
%
% Inputs:
%   accessions - cell array (N x 1) of UniProt accessions (strings). May contain ''.
% Name-value options:
%   'PauseBetweenRequests' (default 0.25) - seconds between HTTP requests
%   'Timeout'              (default 30)   - weboptions timeout (s)
%
% Outputs:
%   foundIds  - Mx1 cell of accession strings that have both mass and sequence
%   MWs       - Mx1 numeric vector of mass values
%   seqs      - Mx1 cell of sequences (strings)
%   foundIdx  - Mx1 numeric indices pointing to positions in the input 'accessions'
%
% Only entries where both mass (numeric) AND sequence (non-empty) were retrieved
% are returned. Order in outputs follows the order of scanning input accessions.

% parse options
p = inputParser;
addParameter(p,'PauseBetweenRequests',0.25,@(x)isnumeric(x) && isscalar(x));
addParameter(p,'Timeout',30,@(x)isnumeric(x) && isscalar(x));
parse(p,varargin{:});
pauseSec = p.Results.PauseBetweenRequests;
timeoutSec = p.Results.Timeout;

% normalize input
if ischar(accessions)
    accessions = {accessions};
end
accessions = accessions(:); % column

n = numel(accessions);
foundIds = {};
MWs = [];
seqs = {};
foundIdx = [];

% weboptions
optsBinary = weboptions('Timeout', timeoutSec, 'ContentType', 'binary'); % download raw bytes

for i = 1:n
    acc = strtrim(accessions{i});
    if isempty(acc)
        continue;
    end

    url = sprintf(['https://rest.uniprot.org/uniprotkb/stream?query=accession:%s' ...
                   '&fields=accession,mass,sequence&format=tsv'], urlencode(acc));

    % --- download to temp file (binary) ---
    tmpFile = [tempname '.bin'];
    try
        websave(tmpFile, url, optsBinary);
    catch ME
        warning('Request failed for %s: %s', acc, ME.message);
        pause(pauseSec);
        if exist(tmpFile,'file'); delete(tmpFile); end
        continue;
    end

    % --- inspect first bytes to detect gzip ---
    fid = fopen(tmpFile, 'rb');
    if fid == -1
        warning('Could not open temp file for %s', acc);
        if exist(tmpFile,'file'); delete(tmpFile); end
        pause(pauseSec);
        continue;
    end
    header = fread(fid, 2, 'uint8')';
    fclose(fid);

    txt = '';
    cleanupFiles = {tmpFile}; % keep track to clean
    try
        if numel(header) >= 2 && header(1) == 31 && header(2) == 139
            % gzip detected -> rename to .gz and gunzip
            gzFile = [tmpFile '.gz'];
            movefile(tmpFile, gzFile);
            cleanupFiles{end} = gzFile;
            try
                outFiles = gunzip(gzFile, tempdir);
            catch
                warning('gunzip failed for %s', gzFile);
                % cleanup and continue
                for f = cleanupFiles, if exist(f{1},'file'), delete(f{1}); end, end
                pause(pauseSec);
                continue;
            end
            if isempty(outFiles)
                warning('gunzip produced no output for %s', gzFile);
                for f = cleanupFiles, if exist(f{1},'file'), delete(f{1}); end, end
                pause(pauseSec);
                continue;
            end
            dataFile = outFiles{1};
            % read as UTF-8 text
            fid2 = fopen(dataFile, 'r', 'n', 'UTF-8');
            if fid2 == -1
                warning('Could not open extracted file %s', dataFile);
            else
                txt = fread(fid2, '*char')';
                fclose(fid2);
            end
            cleanupFiles{end+1} = dataFile;
        else
            % not gzip -> read tmpFile as UTF-8 text
            fid2 = fopen(tmpFile, 'r', 'n', 'UTF-8');
            if fid2 == -1
                warning('Could not open downloaded file %s', tmpFile);
                for f = cleanupFiles, if exist(f{1},'file'), delete(f{1}); end, end
                pause(pauseSec);
                continue;
            end
            txt = fread(fid2, '*char')';
            fclose(fid2);
        end

        % basic check: if returned content looks like HTML or empty, skip
        if isempty(txt)
            warning('Empty response for %s', acc);
            for f = cleanupFiles, if exist(f{1},'file'), delete(f{1}); end, end
            pause(pauseSec);
            continue;
        end
        headerSnippet = lower(txt(1:min(200,end)));
        if contains(headerSnippet, '<html') || contains(headerSnippet, '<!doctype') || contains(headerSnippet, 'error')
            warning('Non-TSV / HTML returned for %s; skipping.', acc);
            for f = cleanupFiles, if exist(f{1},'file'), delete(f{1}); end, end
            pause(pauseSec);
            continue;
        end

        % parse response: split lines and parse last non-empty line
        lines = splitlines(strtrim(txt));
        lines = lines(~cellfun(@(s) isempty(strtrim(s)), lines));
        if numel(lines) < 2
            % no data line found
            for f = cleanupFiles, if exist(f{1},'file'), delete(f{1}); end, end
            pause(pauseSec);
            continue;
        end

        dataLine = lines{end};
        cols = strsplit(dataLine, '\t');
        if numel(cols) < 3
            for f = cleanupFiles, if exist(f{1},'file'), delete(f{1}); end, end
            pause(pauseSec);
            continue;
        end

        accOut   = strtrim(cols{1});
        massStr  = strtrim(cols{2});
        seqStr   = strtrim(cols{3});

        massVal = str2double(massStr);
        if ~isnan(massVal) && ~isempty(seqStr)
            foundIds{end+1,1} = accOut;
            MWs(end+1,1) = massVal;
            seqs{end+1,1} = seqStr;
            foundIdx(end+1,1) = i;
        end

    catch ME
        warning('Error processing response for %s: %s', acc, ME.message);
    end

    % cleanup temp files
    for f = cleanupFiles, if exist(f{1},'file'), delete(f{1}); end, end

    pause(pauseSec); % be polite to API
end

% ensure outputs are columns
foundIds = reshape(foundIds, [], 1);
seqs = reshape(seqs, [], 1);
end
