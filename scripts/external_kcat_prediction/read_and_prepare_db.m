function T = read_and_prepare_db(csvPath)
% Read experimental kcat DB CSV (BRENDA/SABIO merged flavor) and keep only
% columns needed for matching + value. If a substrate-name column exists,
% normalize it to Substrate_norm. If 'ec' exists, keep it as string.

    if ~isfile(csvPath)
        error('BuildKcatMatches:dbMissing','DB file not found: %s', csvPath);
    end

    colsAsString = {'uniprot','InChIKey','MetaNetXID','sequence'};
    keepCols     = {'uniprot','InChIKey','MetaNetXID','sequence','value'};

    nameCandidates = {'Substrate','Substrate_name','substrate','Metabolite','metabolite','name','Name'};

    opts = detectImportOptions(csvPath, 'VariableNamingRule','preserve', 'TextType','string');

    % Optional 'ec'
    if any(strcmp(opts.VariableNames, 'ec'))
        keepCols{end+1} = 'ec';
        colsAsString{end+1} = 'ec';
    end

    % Optional substrate-like
    nameCol = intersect(nameCandidates, opts.VariableNames, 'stable');
    if ~isempty(nameCol)
        keepCols{end+1} = nameCol{1};
    end

    missing = setdiff({'uniprot','InChIKey','MetaNetXID','sequence','value'}, opts.VariableNames);
    if ~isempty(missing)
        error('read_and_prepare_db:MissingColumns','Missing required columns: %s', strjoin(missing, ', '));
    end

    presentStringCols = intersect([colsAsString, nameCol], opts.VariableNames);
    if ~isempty(presentStringCols)
        opts = setvartype(opts, presentStringCols, 'string');
        opts = setvaropts(opts, presentStringCols, 'WhitespaceRule','preserve');
    end

    opts = setvartype(opts, 'value', 'double');
    opts.SelectedVariableNames = keepCols;

    T = readtable(csvPath, opts);

    % Normalize strings
    for i = 1:numel(colsAsString)
        v = colsAsString{i};
        if ismember(v, T.Properties.VariableNames)
            if ~isstring(T.(v)), T.(v) = string(T.(v)); end
            m = ismissing(T.(v)); if any(m), T.(v)(m) = ""; end
        end
    end

    % value_log10
    T.value_log10 = log10(T.value);

    % Substrate_norm
    if ~isempty(nameCol)
        if ~ismember('Substrate', T.Properties.VariableNames)
            T = safe_rename(T, nameCol{1}, 'Substrate');
        end
    else
        T.Substrate = strings(height(T),1);
    end
    T.Substrate_norm = normalize_substrate_name(T.Substrate);
end

function sN = normalize_substrate_name(s)
% Normalize metabolite names for robust equality:
% - lowercase, trim
% - strip trailing bracketed notes
% - unify separators to single space
% - collapse multiple spaces
    if ~isstring(s), s = string(s); end
    sN = lower(strtrim(s));
    sN = regexprep(sN, '\s*(\[[^\]]*\]|\([^\)]*\))\s*$', '', 'once');
    sN = regexprep(sN, '[_\-\,;]+', ' ');
    sN = regexprep(sN, '\s+', ' ');
    sN = strtrim(sN);
end

function T = safe_rename(T, oldName, newName)
    if ismember(oldName, T.Properties.VariableNames) && ~strcmp(oldName, newName)
        try
            T = renamevars(T, oldName, newName); % R2020a+
        catch
            T.(newName) = T.(oldName);
            T.(oldName) = [];
        end
    end
end