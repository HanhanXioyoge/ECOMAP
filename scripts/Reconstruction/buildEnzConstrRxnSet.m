function [rxnWG_kept, report] = buildEnzConstrRxnSet(model, opts)
% buildEnzConstrRxnSet (ECOMAP)
% Build the enzyme-limited reaction set while excluding boundary/transport-like
% reactions using a deterministic rule order. No EC-based dropping is applied.
%
% Removal order (applied in this exact sequence):
%   R0_gene   : drop reactions without gene association (sum(rxnGeneMat,2) == 0)
%   R2_exch   : drop exchange/boundary reactions
%   R3_pure   : drop pure same-entity compartment-switch transport
%   R4_proton : drop proton-coupled translocation ONLY IF:
%                 (i) H+ appears on both sides and crosses compartments, AND
%                (ii) at least one NON-H metabolite also crosses compartments
%   R5_stoich : drop stoichiometric identity (LHS & RHS identical after removing [comp])
%   R6_name   : single lenient name heuristic (merged keywords; case-insensitive),
%               applied only when the reaction has NO EC code (by default)
%
% Report (clear & audit-friendly):
%   report.summary         : struct with total kept/dropped + per-rule counts
%   report.rules           : table listing rule key, description, enabled flag, dropped count
%   report.flags           : table (one row per reaction) with:
%                              rxnID, rxnName, eccodeFirst, ecTopClassFirst, ecHasAny7,
%                              R0_gene, R2_exch, R3_pure, R4_proton, R5_stoich, R6_name, keep
%   report.drop.idx        : struct of index vectors for each rule and 'all'
%   report.drop.tables     : struct of tables listing dropped reactions per rule and 'all'
%
% Options (defaults shown):
%   opts.useNameHeuristic        = true
%   opts.nameHeuristicRequireNoEC= true     % if true, R6 never drops reactions that have EC
%   opts.useProtonCoupling       = true
%   opts.useStoichFallback       = true
%   opts.keywords                = {...}    % merged block-list words (substring match)
%   opts.keepKeywords            = {'pts','phosphotransferase'}  % allow-list
%   opts.extraBlacklist          = {}
%   opts.extraWhitelist          = {}
%   opts.protonRegex             = '^(h|h\+)\[[a-z]\]$'
%   opts.verbose                 = false
%
% Origin: adapted from GECKO-style pre-processing; rewritten and hardened for ECOMAP.

    % ---------- normalize & defaults ----------
    if nargin < 2 || isempty(opts), opts = struct(); end
    if ~isstruct(opts), opts = struct(); end
    if ~isscalar(opts), opts = opts(1); end
    if ~isfield(opts,'verbose') || isempty(opts.verbose), opts.verbose = false; end

    def.useNameHeuristic         = true;
    def.nameHeuristicRequireNoEC = true;   % EC-coded reactions protected from name-based drop
    def.useProtonCoupling        = true;
    def.useStoichFallback        = true;
    def.keywords = {'transport','transporter','translocase','symport','antiport', ...
                    'channel','permease','pump','export','import', ...
                    'uptake','secretion','efflux','influx','shuttle'};
    def.keepKeywords   = {'pts','phosphotransferase'};
    def.extraBlacklist = {};
    def.extraWhitelist = {};
    def.protonRegex    = '^(h|h\+)\[[a-z]\]$';

    fn = fieldnames(def);
    for i = 1:numel(fn)
        f = fn{i};
        if ~isfield(opts,f) || isempty(opts.(f)), opts.(f) = def.(f); end
    end

    % robust normalizer → column string arrays
    toStrVec = @(x) localToStrVec(x);
    opts.keywords       = toStrVec(getFieldOr(opts,'keywords',def.keywords));
    opts.keepKeywords   = toStrVec(getFieldOr(opts,'keepKeywords',def.keepKeywords));
    opts.extraBlacklist = toStrVec(getFieldOr(opts,'extraBlacklist',def.extraBlacklist));
    opts.extraWhitelist = toStrVec(getFieldOr(opts,'extraWhitelist',def.extraWhitelist));

    % sanity checks
    need = {'rxns','rxnNames','S','mets','metNames','rxnGeneMat','lb','ub'};
    for i = 1:numel(need)
        if ~isfield(model, need{i})
            error('model.%s is required.', need{i});
        end
    end

    nRxn     = numel(model.rxns);
    keepMask = true(nRxn,1);

    % ---------- EC diagnostics (for report only; never used to drop) ----------
    [ecFirst, ecTopFirst, ecAny7] = collectRxnECs_(model);
    hasEC = (strlength(strtrim(ecFirst)) > 0) & ~strcmpi(ecFirst, "NA");

    % ---------- R0) gene-associated only ----------
    R0_gene = ~(sum(model.rxnGeneMat, 2) > 0);
    keepMask = keepMask & ~R0_gene;

    % ---------- R2) exchange & R3) pure transport ----------
    [~, exchIdx, ~, pureIdx] = getBoundaryAndTransport(model, 'all');
    R2_exch = false(nRxn,1); R2_exch(exchIdx) = true;
    R3_pure = false(nRxn,1); R3_pure(pureIdx) = true;
    keepMask = keepMask & ~R2_exch & ~R3_pure;

    % helpers
    rxnNamesStr  = lower(string(model.rxnNames(:)));
    rxnText      = rxnNamesStr + " " + lower(string(model.rxns(:))); % match both names & IDs
    S            = model.S;
    metsStr      = string(model.mets(:));
    metNamesStr  = string(model.metNames(:));
    baseMet      = regexprep(metsStr, '\[.*\]$', ''); % strip [comp]
    isProton     = regexpi_all(metsStr, opts.protonRegex) | strcmpi(metNamesStr, 'H+');

    % ---------- R4) proton-coupled translocation ----------
    R4_proton = false(nRxn,1);
    if opts.useProtonCoupling
        for j = 1:nRxn
            if ~keepMask(j), continue; end
            nz  = find(S(:,j) ~= 0); if numel(nz) < 2, continue; end
            neg = nz(S(nz,j) < 0);   pos = nz(S(nz,j) > 0);

            idxHneg = neg(isProton(neg));
            idxHpos = pos(isProton(pos));
            if isempty(idxHneg) || isempty(idxHpos), continue; end

            % (1) H+ must cross compartments
            protonCross = false;
            for a = idxHneg.'
                ca = getCompStr(model, a);
                for b = idxHpos.'
                    cb = getCompStr(model, b);
                    if ~strcmp(ca, cb), protonCross = true; break; end
                end
                if protonCross, break; end
            end
            if ~protonCross, continue; end

            % (2) at least one non-H metabolite must also cross
            nonHCross = false;
            for a = neg(~isProton(neg)).'
                ca = getCompStr(model, a);
                for b = pos(~isProton(pos)).'
                    cb = getCompStr(model, b);
                    if ~strcmp(ca, cb), nonHCross = true; break; end
                end
                if nonHCross, break; end
            end

            R4_proton(j) = protonCross && nonHCross;
        end
        keepMask = keepMask & ~R4_proton;
    end

    % ---------- R5) stoichiometric identity fallback ----------
    R5_stoich = false(nRxn,1);
    if opts.useStoichFallback
        for j = 1:nRxn
            if ~keepMask(j), continue; end
            nz = find(S(:,j) ~= 0);
            if isempty(nz), continue; end
            neg = nz(S(nz,j) < 0); pos = nz(S(nz,j) > 0);
            if isempty(neg) || isempty(pos), continue; end
            if isequal(sort(baseMet(neg)), sort(baseMet(pos)))
                R5_stoich(j) = true;
            end
        end
        keepMask = keepMask & ~R5_stoich;
    end

    % ---------- R6) single lenient name heuristic ----------
    R6_name = false(nRxn,1);
    if opts.useNameHeuristic
        blockMask = false(nRxn,1);
        allowMask = false(nRxn,1);

        % merged blacklist
        for kw = [opts.keywords; opts.extraBlacklist].'
            blockMask = blockMask | contains(rxnText, lower(kw));
        end
        % whitelist
        for kw = [opts.keepKeywords; opts.extraWhitelist].'
            allowMask = allowMask | contains(rxnText, lower(kw));
        end
        blockMask(allowMask) = false; % whitelist wins

        % optional EC protection
        if opts.nameHeuristicRequireNoEC
            blockMask = blockMask & ~hasEC;
        end

        R6_name = blockMask;
        keepMask = keepMask & ~R6_name;
    end

    % ---------- outputs ----------
    rxnWG_kept = find(keepMask);

    % ============ REPORT ============
    report.flags = table( ...
        model.rxns(:), model.rxnNames(:), ...
        ecFirst, ecTopFirst, ecAny7, ...
        R0_gene, R2_exch, R3_pure, R4_proton, R5_stoich, R6_name, ...
        keepMask, ...
        'VariableNames', {'rxnID','rxnName','eccodeFirst','ecTopClassFirst','ecHasAny7', ...
                          'R0_gene','R2_exch','R3_pure','R4_proton','R5_stoich','R6_name','keep'} ...
    );

    report.drop.idx = struct();
    report.drop.idx.R0_gene   = find(R0_gene);
    report.drop.idx.R2_exch   = find(R2_exch);
    report.drop.idx.R3_pure   = find(R3_pure);
    report.drop.idx.R4_proton = find(R4_proton);
    report.drop.idx.R5_stoich = find(R5_stoich);
    report.drop.idx.R6_name   = find(R6_name);
    report.drop.idx.all       = find(~keepMask);

    mkTable = @(idx) table( ...
        model.rxns(idx), model.rxnNames(idx), ...
        ecFirst(idx), ecTopFirst(idx), ecAny7(idx), ...
        'VariableNames', {'rxnID','rxnName','eccodeFirst','ecTopClassFirst','ecHasAny7'} ...
    );

    report.drop.tables = struct();
    report.drop.tables.R0_gene   = mkTable(report.drop.idx.R0_gene);
    report.drop.tables.R2_exch   = mkTable(report.drop.idx.R2_exch);
    report.drop.tables.R3_pure   = mkTable(report.drop.idx.R3_pure);
    report.drop.tables.R4_proton = mkTable(report.drop.idx.R4_proton);
    report.drop.tables.R5_stoich = mkTable(report.drop.idx.R5_stoich);
    report.drop.tables.R6_name   = mkTable(report.drop.idx.R6_name);
    report.drop.tables.all       = mkTable(report.drop.idx.all);

    report.summary = struct();
    report.summary.totalRxn  = full(double(nRxn));
    report.summary.totalKeep = full(double(sum(keepMask)));
    report.summary.totalDrop = full(double(sum(~keepMask)));

    ruleKeys   = {'R0_gene','R2_exch','R3_pure','R4_proton','R5_stoich','R6_name'};
    ruleDesc   = { ...
        'Not gene-associated', ...
        'Exchange/boundary', ...
        'Pure same-entity compartment switch', ...
        'Proton-coupled (H+ crosses AND at least one non-H crosses)', ...
        'Stoichiometric identity (LHS = RHS)', ...
        'Name heuristic (keywords/whitelist)' ...
    };
    ruleEnabled = [true, true, true, opts.useProtonCoupling, opts.useStoichFallback, opts.useNameHeuristic];
    droppedCount = cellfun(@(k) numel(report.drop.idx.(k)), ruleKeys);

    report.summary.byRule = cell2struct(num2cell(droppedCount), ruleKeys, 2);
    report.rules = table( ...
        ruleKeys(:), ruleDesc(:), logical(ruleEnabled(:)), droppedCount(:), ...
        'VariableNames', {'RuleKey','Description','Enabled','DroppedCount'} ...
    );

    if opts.verbose
        fprintf('[buildEnzConstrRxnSet] kept=%d, dropped=%d of %d (nameHeuristicRequireNoEC=%s)\n', ...
            report.summary.totalKeep, report.summary.totalDrop, nRxn, string(opts.nameHeuristicRequireNoEC));
    end
end

% ---------- helpers ----------
function tf = regexpi_all(strCell, pattern)
    if isstring(strCell), strCell = cellstr(strCell); end
    tf = false(numel(strCell),1);
    for i=1:numel(strCell)
        tf(i) = ~isempty(regexpi(strCell{i}, pattern, 'once'));
    end
end

function val = getFieldOr(S, name, fallback)
    if isfield(S, name) && ~isempty(S.(name))
        val = S.(name);
    else
        val = fallback;
    end
end

function s = localToStrVec(x)
    if isempty(x)
        s = strings(0,1);
    elseif isstring(x)
        s = x(:);
    elseif iscellstr(x) || (iscell(x) && all(cellfun(@ischar,x)))
        s = string(x(:));
    elseif ischar(x)
        s = string({x});
    else
        s = string(x(:));
    end
end

function [ecFirst, ecTopFirst, ecAny7] = collectRxnECs_(model)
% Only use model.eccodes (auditing only; never used to drop).
% ecFirst    : first EC token per reaction ("" if none)
% ecTopFirst : first EC top class digit (1..7) or NaN
% ecAny7     : true if ANY token is 7.x (otherwise false)
    n = numel(model.rxns);
    ecFirst    = strings(n,1);
    ecTopFirst = nan(n,1);
    ecAny7     = false(n,1);

    if isfield(model,'eccodes') && numel(model.eccodes) == n
        ecRaw = string(model.eccodes(:));
    else
        ecRaw = strings(n,1); % no ECs available
    end

    for i = 1:n
        s = strtrim(ecRaw(i));
        if s == "" || strcmpi(s,"NA"), continue; end
        toks = string(regexp(s, '[^;, ]+', 'match'));
        if isempty(toks), continue; end

        % first token
        ecFirst(i) = toks(1);
        m = regexp(ecFirst(i), '^([0-9])\.', 'tokens','once');
        if ~isempty(m)
            ecTopFirst(i) = str2double(m{1});
        elseif startsWith(ecFirst(i), "7.")
            ecTopFirst(i) = 7;
        end

        % any 7.x (for auditing)
        ecAny7(i) = any(startsWith(toks, "7."));
    end
end

function c = getCompStr(model, mIdx)
    % Return compartment label for a metabolite index (as string).
    if isfield(model,'metComps') && ~isempty(model.metComps)
        c = num2str(model.metComps(mIdx));  % could map to names if needed
    else
        mID = model.mets{mIdx};
        tok = regexp(mID, '\[([^\]]+)\]$', 'tokens', 'once');
        if isempty(tok), c = ''; else, c = tok{1}; end
    end
end
