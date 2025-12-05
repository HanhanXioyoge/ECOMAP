function [kcatList_out, idxKept, idxDropped, stats] = remapKcatListByType(nameMap, typeFrom, typeTo, kcatList)
% REMAPKCATLISTBYTYPE
% Map kcatList.rxns names from one ecModel type to another using a nameMap table
% and drop any rows that do not have a mapping. Aligned fields (genes,
% substrates, kcats) are trimmed at the same indices.
%
% Inputs
%   nameMap   : table produced by your intersection mapper (columns like: basic, integrated, isozyme)
%   typeFrom  : source type name, e.g., 'basic' | 'integrated' | 'isozyme'
%   typeTo    : destination type name (same set as above)
%   kcatList  : struct with at least field 'rxns'. Optional aligned fields:
%               - kcatList.genes
%               - kcatList.substrates
%               - kcatList.kcats
%               These fields must be vectors or matrices whose first dimension
%               (or numel for vectors) equals numel(kcatList.rxns).
%
% Outputs
%   kcatList_out : same struct as input but with rxn names replaced to typeTo
%                  and non-mapped rows removed across aligned fields.
%   idxKept      : indices (in original kcatList) that were kept/mapped
%   idxDropped   : indices that were removed (no mapping)
%   stats        : struct with counts (nIn, nMapped, nDropped)
%
% Notes
% - Column picking is done by sanitized variable names:
%     col = matlab.lang.makeValidName(lower(type))
%   which should match how you built nameMap ('basic','integrated','isozyme').
% - Mapping is exact on the strings in nameMap.(typeFrom) → nameMap.(typeTo).
% - If kcatList.rxns contains duplicates, each occurrence is mapped independently.
% - If kcatList.kcats is a matrix with N rows (N = numel(rxns)), we keep rows by mask.

    % -------- normalize & validate inputs --------
    if nargin < 4
        error('Usage: remapKcatListByType(nameMap, typeFrom, typeTo, kcatList)');
    end
    if ~istable(nameMap)
        error('nameMap must be a table with columns for ecModel types (e.g., basic, integrated, isozyme).');
    end
    if ~isfield(kcatList, 'rxns')
        error('kcatList must contain a field ''rxns''.');
    end

    % Column names to pick (sanitized, lowercased)
    colFrom = char(matlab.lang.makeValidName(lower(string(typeFrom))));
    colTo   = char(matlab.lang.makeValidName(lower(string(typeTo))));

    if ~ismember(colFrom, nameMap.Properties.VariableNames)
        error('nameMap does not contain source type column "%s".', colFrom);
    end
    if ~ismember(colTo, nameMap.Properties.VariableNames)
        error('nameMap does not contain destination type column "%s".', colTo);
    end

    % -------- build the mapping arrays --------
    srcNames = string(nameMap.(colFrom));
    dstNames = string(nameMap.(colTo));

    % Remove rows with empty src or empty dst (safety)
    ok = (strlength(srcNames) > 0) & (strlength(dstNames) > 0);
    srcNames = srcNames(ok);
    dstNames = dstNames(ok);

    % -------- map kcatList.rxns --------
    rxnsIn      = kcatList.rxns;
    wasCell     = iscell(rxnsIn);
    rxnsIn_str  = string(rxnsIn(:));   % column string array
    N           = numel(rxnsIn_str);

    % ismember gives locations into srcNames
    [tf, loc] = ismember(rxnsIn_str, srcNames);
    idxKept    = find(tf);
    idxDropped = find(~tf);

    % mapped names for kept rows
    newNames = dstNames(loc(tf));

    % -------- assemble output kcatList --------
    kcatList_out = kcatList;  % start with a copy

    % Replace & trim rxns
    if wasCell
        rxnOut = cell(N,1);
        rxnOut(idxKept) = cellstr(newNames);
        rxnOut = rxnOut(idxKept);
    else
        rxnOut = strings(N,1);
        rxnOut(idxKept) = newNames;
        rxnOut = rxnOut(idxKept);
    end
    kcatList_out.source = kcatList.source;
    kcatList_out.rxns = rxnOut;

    % Trim aligned fields if present
    kcatList_out = tryTrimField(kcatList_out, 'genes',      idxKept, N);
    kcatList_out = tryTrimField(kcatList_out, 'substrates', idxKept, N);
    kcatList_out = tryTrimField(kcatList_out, 'kcats',      idxKept, N);

    % -------- stats --------
    stats = struct();
    stats.nIn      = N;
    stats.nMapped  = numel(idxKept);
    stats.nDropped = numel(idxDropped);
end

% ====================== local helpers ======================

function S = tryTrimField(S, fieldName, keepIdx, N)
% Trim S.(fieldName) to keepIdx if it appears aligned with rxns
% Rules:
%   - If the field is a vector with numel == N -> index by keepIdx
%   - If the field is a matrix with size(field,1) == N -> keep rows keepIdx
%   - Otherwise leave it untouched

    if ~isfield(S, fieldName) || isempty(S.(fieldName))
        return;
    end

    val = S.(fieldName);

    % vector case (row or column)
    if isvector(val) && numel(val) == N
        S.(fieldName) = val(keepIdx);
        return;
    end

    % matrix/table timetable: use first dimension if it equals N
    if (ismatrix(val) || istable(val)) && size(val,1) == N
        S.(fieldName) = val(keepIdx, :);
        return;
    end

    % cell array with N rows
    if iscell(val) && size(val,1) == N
        S.(fieldName) = val(keepIdx, :);
        return;
    end

    % If none of the above match, we skip trimming to avoid shape errors.
end
