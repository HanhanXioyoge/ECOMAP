function cleanName = cleanEnzymeName(name, opts)
%CLEANENZYMENAME Whitelist-based stripping of model bookkeeping suffixes.
%
%   cleanName = cleanEnzymeName(name) strips model bookkeeping suffixes
%   from a reaction/enzyme name so the result is suitable for UniProt
%   protein_name: searches.
%
%   Default whitelist (bookkeeping only):
%     - Trailing parenthetical tokens: (No1), (No2), (No3), (arm), ...
%     - Trailing underscore tokens: _REV, _No1, _No2, _No3
%     - Leading arm_ prefix
%
%   NOT stripped (functional definition, must be preserved):
%     - (ATP), (NAD+), (NADP+), (NADPH), (CoA), (decarboxylating), ...
%     - [NADP+], [NAD+], [NADPH]
%
%   cleanName = cleanEnzymeName(name, opts) accepts:
%     opts.StripSuffixes (cellstr) - additional whitelist entries
%
%   See: docs/superpowers/specs/2026-08-09-oko-plus-homolog-interval-pipeline-design.md
%        section 4.1

    if nargin < 2
        opts = struct();
    end

    if ~isfield(opts, 'StripSuffixes')
        opts.StripSuffixes = {};
    end

    % 默认 whitelist：仅模型 bookkeeping 后缀（严格匹配 spec §4.1）
    % 'arm_' 仅作为前导前缀处理，不在此处作为 trailing suffix 出现，
    % 否则会被末尾匹配错误地 strip（如 'something_arm_' -> 'something'）。
    defaultWhitelist = {
        '(No1)', '(No2)', '(No3)', ...
        '(arm)', '(No1)_reverse', '(No2)_reverse', '(No3)_reverse', ...
        '_REV', '_No1', '_No2', '_No3'
    };

    allSuffixes = [defaultWhitelist(:); opts.StripSuffixes(:)];

    cleanName = strtrim(name);

    % 1. 去除前导 arm_（特殊处理）
    if startsWith(cleanName, 'arm_')
        cleanName = cleanName(5:end);
    end

    % 2. 末尾后缀剥离（按最长匹配优先，避免短后缀吃掉长后缀的部分）
    %    spec §4.1：最长匹配优先；先按 suffix 长度降序排序，保证重叠场景下
    %    总是匹配最长的那个（如 opts.StripSuffixes = {'reverse','_No1_reverse'}
    %    对 'X_No1_reverse' 应得到 'X'）。
    if ~isempty(allSuffixes)
        [~, sortIdx] = sort(cellfun(@length, allSuffixes), 'descend');
        allSuffixes = allSuffixes(sortIdx);
    end

    matched = true;
    while matched
        matched = false;
        for i = 1:numel(allSuffixes)
            suffix = allSuffixes{i};
            if endsWith(cleanName, suffix)
                cleanName = strtrim(cleanName(1:end-length(suffix)));
                matched = true;
                break;
            end
        end
    end

    % 3. 合并多余空格
    cleanName = regexprep(cleanName, '\s+', ' ');
end
