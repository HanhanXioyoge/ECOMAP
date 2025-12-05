function [model, noMNX] = addMetMetaNetXID(model, verbosity)
% addMetMetaNetXID
% Populate per-metabolite MetaNetX Chemical IDs (MNXM...) into model.metMetaNetXID
% by parsing model.metMiriams. "Virtual" metabolites whose names start with
% 'prot_' are ignored for both extraction and coverage (name is blanked).
%
% INPUTS
%   model     : struct with
%               - metMiriams : Nx1 cell; each cell is a struct or a cell-array of structs
%                              containing cross-references, with fields .name / .value
%               - metNames   : Nx1 cellstr/strings (used for prot_-filtering and reporting)
%               - (optional) metMetaNetXID : Nx1 cellstr; existing non-empty values are preserved
%   verbosity : 0 = silent, 1 = info (default), 2 = debug
%
% OUTPUTS
%   model     : with model.metMetaNetXID (Nx1 cellstr), each entry either 'MNXM####' or ''
%   noMNX     : list of metabolite names that are still empty AFTER ignoring prot_* entries
%
% BEHAVIOR
%   - Names preprocessing: remove any name matching '^prot_.*' by blanking the string,
%     then strtrim. Blank names are "ignored" for extraction & coverage.
%   - Extraction prioritizes MIRIAM entries whose .name contains 'metanetx.chemical';
%     then falls back to scanning all .value fields for tokens 'MNXM\d+'.
%   - Existing non-empty model.metMetaNetXID entries are kept (fill-missing-only).
%
% NOTE
%   - If metNames is missing, the function cannot identify prot_* entries and will not
%     ignore any metabolites for coverage; a warning is printed at verbosity>=1.

    if nargin < 2 || isempty(verbosity), verbosity = 1; end

    % ----- size and names -----
    if isfield(model,'metNames') && ~isempty(model.metNames)
        names = model.metNames;
        if isstring(names), names = cellstr(names); end
        names = names(:);
    else
        % If names are missing, synthesize placeholders for reporting only
        if verbosity >= 1
            fprintf('[addMetMetaNetXID] Warning: model.metNames not found; cannot apply prot_* filtering.\n');
        end
        if isfield(model,'metMiriams') && ~isempty(model.metMiriams)
            N = numel(model.metMiriams);
        else
            error('addMetMetaNetXID: model.metMiriams or model.metNames is required.');
        end
        names = arrayfun(@(i) sprintf('#%d',i), (1:N).', 'UniformOutput', false);
    end

    N = numel(names);

    % ----- prot_* filtering (ignore virtual metabolites) -----
    % Blank any name that matches '^prot_.*' and trim spaces
    names = regexprep(names, '^prot_.*', '');
    names = cellfun(@strtrim, names, 'UniformOutput', false);
    isIgnored = cellfun(@isempty, names);  % ignored in extraction & coverage

    % ----- ensure target field dimensions -----
    if ~isfield(model,'metMetaNetXID') || isempty(model.metMetaNetXID)
        model.metMetaNetXID = repmat({''}, N, 1);
    else
        v = model.metMetaNetXID;
        if isstring(v), v = cellstr(v); end
        v = v(:);
        if numel(v) ~= N
            vv = repmat({''}, N, 1);
            k = min(N, numel(v));
            vv(1:k) = v(1:k);
            v = vv;
        end
        model.metMetaNetXID = v;
    end

    % ----- main loop: fill missing MNX IDs (skip ignored rows) -----
    filled = 0;
    for i = 1:N
        if isIgnored(i), continue; end                           % skip prot_*
        if ~isempty(orEmpty(model.metMetaNetXID{i})), continue; end % keep existing non-empty

        if ~isfield(model,'metMiriams') || i > numel(model.metMiriams) || isempty(model.metMiriams{i})
            continue;
        end

        miriam = model.metMiriams{i};
        ids = extract_mnx_ids_from_miriam(miriam);  % {'MNXM####', ...}

        if ~isempty(ids)
            model.metMetaNetXID{i} = ids{1};  % take the first, order is stable
            filled = filled + 1;
            if verbosity >= 2
                fprintf('[MNX][%d] %s -> %s\n', i, safeName(names{i}), model.metMetaNetXID{i});
            end
        end
    end

    % ----- coverage (exclude ignored rows) -----
    emptyMask = cellfun(@isempty, model.metMetaNetXID) & ~isIgnored;
    denom = sum(~isIgnored);
    got = sum(~emptyMask & ~isIgnored);
    cov = ifthen(denom==0, 0, got / denom);

    if verbosity >= 1
        fprintf('MetaNetX coverage : %.1f%% (%d/%d)\n', 100*cov, got, denom);
        if verbosity >= 2
            fprintf('Ignored (prot_* or blank) metabolites: %d/%d\n', sum(isIgnored), N);
        end
    end

    % ----- output noMNX list (exclude ignored rows) -----
    noMNX = names(emptyMask);
end

% =================== helpers ===================

function ids = extract_mnx_ids_from_miriam(miriam)
% Extract tokens 'MNXM\d+' from a MIRIAM entry/entries.
% Priority: entries whose .name contains 'metanetx.chemical'; then scan all .value.
%
% Accepts:
%   - struct with fields .name/.value
%   - cell array of such structs
%   - any nested shape commonly seen in models
    ids = {};
    mNames = field_to_cellstr(miriam, 'name');
    mVals  = field_to_cellstr(miriam, 'value');

    % prioritize fields explicitly referencing MetaNetX Chemical
    pri = find(contains(lower(mNames), 'metanetx.chemical'));
    for p = pri(:).'
        tok = regexp(mVals{p}, 'MNXM\d+', 'match');
        if ~isempty(tok), ids = [ids, tok]; end %#ok<AGROW>
    end

    % fallback: scan all values for MNXM#### tokens
    if isempty(ids)
        for j = 1:numel(mVals)
            tok = regexp(mVals{j}, 'MNXM\d+', 'match');
            if ~isempty(tok), ids = [ids, tok]; end %#ok<AGROW>
        end
    end

    if ~isempty(ids), ids = unique(ids, 'stable'); end
end

function arr = field_to_cellstr(s, fname)
% Collect a field from a struct or cell-of-structs into a cellstr array.
    arr = {};
    try
        if iscell(s)
            tmp = {};
            for i = 1:numel(s)
                if isstruct(s{i}) && isfield(s{i}, fname)
                    v = s{i}.(fname);
                    if iscell(v)
                        tmp = [tmp, cellfun(@char, v, 'UniformOutput', false)]; %#ok<AGROW>
                    else
                        tmp = [tmp, {char(v)}]; %#ok<AGROW>
                    end
                end
            end
            arr = tmp;
        elseif isstruct(s)
            if isfield(s, fname)
                v = s.(fname);
                if iscell(v)
                    arr = cellfun(@char, v, 'UniformOutput', false);
                else
                    arr = {char(v)};
                end
            end
        end
    catch
        arr = {};
    end
end

function s = orEmpty(x)
% Convert various "empty-like" inputs to '' (char), otherwise char(string(x)).
    try
        if isempty(x), s = ''; return; end
        if iscell(x)
            if isempty(x{1}), s = ''; else, s = orEmpty(x{1}); end
            return;
        end
        if isstring(x)
            x = x(:);
            mask = ~ismissing(x) & (strlength(x) > 0);
            if any(mask), s = char(x(find(mask,1,'first'))); else, s = ''; end
            return;
        end
        if ischar(x), s = x; return; end
        s = char(string(x));
        if all(strlength(string(s)) == 0), s = ''; end
    catch
        s = '';
    end
end

function out = safeName(x)
% Shorten long names in logs for readability.
    out = orEmpty(x);
    if numel(out) > 60, out = [out(1:60) '...']; end
end

function v = ifthen(cond, a, b)
% Ternary helper.
    if cond, v = a; else, v = b; end
end