function report = test_chebi_by_id(chebiIds, varargin)
% test_chebi_by_id.m
% Query ChEBI 2.0 "public compound" endpoint by CHEBI IDs
%   e.g. https://www.ebi.ac.uk/chebi/backend/api/public/compound/CHEBI:16750/
% Extract Canonical SMILES and InChIKey, with retries and robust JSON parsing.
%
% Usage:
%   R = test_chebi_by_id("CHEBI:16750");
%   R = test_chebi_by_id(["CHEBI:16750","CHEBI:17234"], 'Timeout',15, 'Retries',3);

    if nargin < 1 || isempty(chebiIds)
        chebiIds = "CHEBI:16750";   % default example
    end
    chebiIds = string(chebiIds);

    p = inputParser;
    addParameter(p,'Timeout', 10, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p,'Retries', 2, @(x) isnumeric(x) && isscalar(x) && x>=0);
    addParameter(p,'Verbose', true, @(x) islogical(x) || x==0 || x==1);
    addParameter(p,'Sleep', 0.2, @(x) isnumeric(x) && isscalar(x) && x>=0);
    parse(p, varargin{:});

    Timeout = p.Results.Timeout;
    Retries = p.Results.Retries;
    Verbose = logical(p.Results.Verbose);
    PauseS  = p.Results.Sleep;

    UA = 'MATLAB-ChEBI-ByID/2.0';
    wopt_text = weboptions('Timeout',Timeout, ...
                           'HeaderFields',{'Accept','application/json';'User-Agent',UA}, ...
                           'ContentType','text');

    base = "https://www.ebi.ac.uk/chebi/backend/api/public/compound/";

    if Verbose
        fprintf('[ChEBI public/compound] IDs: %s\n', strjoin(chebiIds, ', '));
        fprintf('[ChEBI public/compound] Timeout=%.1fs Retries=%d\n', Timeout, Retries);
    end

    report = repmat(struct( ...
        'ChebiId',"",'Status',0,'Success',false,'DurationSec',NaN, ...
        'CanonicalSMILES',"",'InChIKey',"",'ErrorMessage',"",'UrlTried',""), ...
        numel(chebiIds), 1);

    for i = 1:numel(chebiIds)
        cid = strtrim(chebiIds(i));
        t0 = tic;
        if Verbose, fprintf('--- [%d/%d] %s\n', i, numel(chebiIds), cid); end

        url = base + encode_id_for_path(cid) + "/";   % 按你给的形式保留末尾斜杠
        ok = false; smiles = ""; ik = ""; lastErr = ""; status = 0;

        try
            txt = tryReadText(url, wopt_text, Retries);
            status = 200;
            J = safeJsonDecode(txt);
            [smiles, ik] = extractSmilesInChIKey(J, txt);
            if smiles~="" || ik~="", ok = true; else, lastErr = "No SMILES/InChIKey in response"; end
        catch ME
            lastErr = ME.message;
        end

        report(i).ChebiId         = cid;
        report(i).Status          = status;
        report(i).Success         = ok;
        report(i).DurationSec     = toc(t0);
        report(i).CanonicalSMILES = string(smiles);
        report(i).InChIKey        = string(ik);
        report(i).ErrorMessage    = string(lastErr);
        report(i).UrlTried        = string(url);

        if Verbose
            if ok
                fprintf('  ✓ %s\n', cid);
                if smiles~="", fprintf('    SMILES   : %s\n', truncate(smiles)); end
                if ik~="",     fprintf('    InChIKey : %s\n', ik); end
            else
                fprintf('  ! %s failed: %s\n', cid, lastErr);
            end
        end

        pause(PauseS);
    end

    if Verbose
        nOK = nnz([report.Success]);
        fprintf('\n[Summary] %d/%d successful.\n', nOK, numel(report));
    end
end

% -------- helpers --------
function txt = tryReadText(url, wopt, retries)
    lastErr = [];
    for a=0:retries
        try
            txt = webread(url, wopt);  % return text
            if ~(ischar(txt) || isstring(txt)), txt = jsonencode(txt); end
            return
        catch ME
            lastErr = ME; pause(0.4);
        end
    end
    if isempty(lastErr), error('GET %s failed (unknown)', url);
    else, rethrow(lastErr);
    end
end

function J = safeJsonDecode(txt)
    J = [];
    try
        J = jsondecode(char(txt));
    catch
        % keep empty; regex fallback will handle
    end
end

function [smiles, inchikey] = extractSmilesInChIKey(J, rawTxt)
    smiles = ""; inchikey = "";
    % 先遍历结构体找关键字段（不依赖具体 schema）
    if ~isempty(J) && isstruct(J)
        [smiles, inchikey] = walkStructForKeys(J);
        if smiles~="" || inchikey~="", return; end
    end
    % 再用正则兜底（兼容大小写和命名变体）
    s = char(string(rawTxt));
    m1 = regexp(s,'"(?:canonical_?smiles|smiles)"\s*:\s*"([^"]+)"','tokens','once');
    if ~isempty(m1), smiles = string(m1{1}); end
    m2 = regexp(s,'"(?:standard_?inchi_?key|inchi_?key)"\s*:\s*"([^"]+)"','tokens','once');
    if ~isempty(m2), inchikey = string(m2{1}); end
end

function [sm, ik] = walkStructForKeys(S)
    sm = ""; ik = "";
    if ~isstruct(S), return; end
    fns = fieldnames(S);
    for i=1:numel(fns)
        fn = fns{i};
        val = S.(fn);
        lfn = lower(fn);
        if ischar(val) || isstring(val)
            if contains(lfn,'smiles') && sm=="", sm = string(val); end
            if contains(lfn,'inchi') && contains(lfn,'key') && ik=="", ik = string(val); end
        elseif isstruct(val)
            [sm2, ik2] = walkStructForKeys(val);
            if sm=="", sm = sm2; end
            if ik=="", ik = ik2; end
        elseif iscell(val)
            for k=1:numel(val)
                [sm2, ik2] = walkStructForKeys(val{k});
                if sm=="", sm = sm2; end
                if ik=="", ik = ik2; end
            end
        end
        if sm~="" && ik~="", return; end
    end
end

function s = encode_id_for_path(idstr)
    % 对路径段做编码（特别是冒号）
    s = string(idstr);
    s = strrep(s, ":", "%3A");
    s = regexprep(s, ' ', '%20');
end

function s = truncate(x)
    x = string(x);
    if strlength(x) > 80, s = extractBefore(x,77) + "..."; else, s = x; end
end
