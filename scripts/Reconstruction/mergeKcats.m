function mergedKcatList = mergeKcats(kcatList_complete, kcatList_prediction, kcatList_fuzzy, topOriginLimit, bottomOriginLimit, wildcardLimit)
% mergeKcats (ECOMAP)
%   Merge three kcat sources into one prioritized list:
%     1) COMPLETE (e.g., curated/experimental matches produced by completeKcatMatch)
%     2) BRENDA FUZZY (exact EC & good origin first; then lower-priority fuzzy)
%     3) PREDICTION (DLKcat/UniKP/CatPred-style)
%
%   Inputs
%     kcatList_complete    : struct with fields rxns, genes, substrates, kcats, source (e.g., 'CompleteMatch')
%     kcatList_prediction  : struct with fields rxns, genes, substrates, kcats, source (e.g., 'DLKcat'/'UniKP'/'CatPred')
%     kcatList_fuzzy       : struct with fields rxns, substrates, kcats, eccodes, wildcardLvl, origin, source ('brenda')
%     topOriginLimit       : (opt, default 6) exact-EC fuzzy with origin ≤ topOriginLimit are high-priority
%     bottomOriginLimit    : (opt, default 6) lower-priority fuzzy need origin ≤ bottomOriginLimit
%     wildcardLimit        : (opt, default 3) max EC wildcard level allowed for lower-priority fuzzy
%
%   Output
%     mergedKcatList       : struct with fields
%         .source       = 'Merged complete + fuzzy + prediction'
%         .kcatSource   : per-row provenance string
%         .rxns, .genes, .substrates, .kcats
%         .eccodes, .wildcardLvl, .origin (NaN where not applicable)
%
%   Priority rules (per reaction):
%     A) Use COMPLETE if available (assumed highest confidence).
%     B) Else, use FUZZY with exact EC (wildcardLvl==0) and origin ≤ topOriginLimit.
%     C) Else, use PREDICTION.
%     D) Else, use remaining FUZZY:
%           exact EC with topOriginLimit < origin ≤ bottomOriginLimit, OR
%           wildcard EC with wildcardLvl ≤ wildcardLimit and origin ≤ bottomOriginLimit.
%
%   Note
%     This function is adapted from GECKO’s merging logic and generalized
%     for ECOMAP. Field names are kept consistent with ECOMAP utilities.

    % -------------------- defaults & checks --------------------
    if nargin < 6, wildcardLimit    = 3; end
    if nargin < 5, bottomOriginLimit = 6; end
    if nargin < 4, topOriginLimit    = 6; end

    if (topOriginLimit < 1) || (topOriginLimit > 6)
        error('topOriginLimit should be between 1 and 6.');
    end
    if (bottomOriginLimit < 1) || (bottomOriginLimit > 6)
        error('bottomOriginLimit should be between 1 and 6.');
    end
    if (wildcardLimit < 0) || (wildcardLimit > 3)
        error('wildcardLimit should be between 0 and 3.');
    end

    % Helper to safely fetch field or fill with NaNs/cells
    function v = getfield_or(k, name, n, fill)
        if isfield(k, name)
            v = k.(name);
        else
            if ischar(fill), v = repmat({fill}, n, 1);
            elseif isnan(fill), v = nan(n,1);
            else, v = repmat(fill, n, 1);
            end
        end
    end

    % -------------------- normalize FUZZY meta --------------------
    wc = getfield_or(kcatList_fuzzy, 'wildcardLvl', numel(kcatList_fuzzy.rxns), NaN);
    orgn = getfield_or(kcatList_fuzzy, 'origin',     numel(kcatList_fuzzy.rxns), NaN);

    wc_norm   = wc;   wc_norm(isnan(wc_norm))   = 1000; % treat missing as very large
    orgn_norm = orgn; orgn_norm(isnan(orgn_norm)) = 1000;

    % High-priority fuzzy (exact EC & good origin)
    prioFuzzy1 = (wc_norm == 0) & (orgn_norm <= topOriginLimit);
    rxnsFuzzy1 = unique(kcatList_fuzzy.rxns(prioFuzzy1));

    % Lower-priority fuzzy mask (after prediction fallback)
    prioFuzzyLow = ((wc_norm == 0) & (orgn_norm > topOriginLimit) & (orgn_norm <= bottomOriginLimit)) | ...
                   ((wc_norm  > 0) & (wc_norm <= wildcardLimit)   & (orgn_norm <= bottomOriginLimit));

    % -------------------- reaction-level sets --------------------
    rxnsComplete   = unique(getfield_or(kcatList_complete,  'rxns', 0, { }));
    rxnsPredict    = unique(getfield_or(kcatList_prediction,'rxns', 0, { }));
    rxnsFuzzyAll   = unique(getfield_or(kcatList_fuzzy,     'rxns', 0, { }));

    % -------------------- build merged list --------------------
    mergedKcatList          = struct();
    mergedKcatList.source   = 'Merged complete + fuzzy + prediction';

    % Collect rows in priority order:
    % A) COMPLETE
    selA = true(numel(rxnsComplete),1); %#ok<NASGU> (not used directly, just semantic)

    % Find fuzzy exact-good that are NOT already covered by complete
    notInComplete = ~ismember(kcatList_fuzzy.rxns, rxnsComplete);
    selB = prioFuzzy1 & notInComplete;

    % For prediction, exclude reactions already in A or B
    rxnsCoveredAB = unique([rxnsComplete; kcatList_fuzzy.rxns(selB)]);
    selC = ~ismember(kcatList_prediction.rxns, rxnsCoveredAB);

    % Remaining fuzzy (low-priority), exclude anything covered by A or C (and B, implicitly)
    rxnsCoveredABC = unique([rxnsCoveredAB; kcatList_prediction.rxns(selC)]);
    selD = prioFuzzyLow & ~ismember(kcatList_fuzzy.rxns, rxnsCoveredABC);

    % Assemble fields with proper NA fill for non-fuzzy rows
    % A) COMPLETE
    A_rxns       = getfield_or(kcatList_complete,  'rxns',       0, { });
    A_genes      = getfield_or(kcatList_complete,  'genes',      numel(A_rxns), '');
    A_subs       = getfield_or(kcatList_complete,  'substrates', numel(A_rxns), '');
    A_kcats      = getfield_or(kcatList_complete,  'kcats',      numel(A_rxns), NaN);
    A_sourceRow  = getfield_or(kcatList_complete,  'source',     1,   'CompleteMatch');
    A_source     = repmat({A_sourceRow}, numel(A_rxns), 1);
    A_eccodes    = repmat({[]}, numel(A_rxns), 1);
    A_wc         = nan(numel(A_rxns),1);
    A_orgn       = nan(numel(A_rxns),1);

    % B) FUZZY (prio1)
    B_rxns       = kcatList_fuzzy.rxns(selB);
    B_genes      = repmat({''}, sum(selB), 1);
    B_subs       = kcatList_fuzzy.substrates(selB);
    B_kcats      = kcatList_fuzzy.kcats(selB);
    B_source     = repmat({kcatList_fuzzy.source}, sum(selB), 1);
    B_eccodes    = getfield_or(kcatList_fuzzy, 'eccodes', numel(kcatList_fuzzy.rxns), '');
    B_eccodes    = B_eccodes(selB);
    B_wc         = wc(selB);
    B_orgn       = orgn(selB);

    % C) PREDICTION
    C_rxns       = kcatList_prediction.rxns(selC);
    C_genes      = getfield_or(kcatList_prediction, 'genes',      numel(kcatList_prediction.rxns), '');
    C_genes      = C_genes(selC);
    C_subs       = getfield_or(kcatList_prediction, 'substrates', numel(kcatList_prediction.rxns), '');
    C_subs       = C_subs(selC);
    C_kcats      = kcatList_prediction.kcats(selC);
    C_source     = repmat({kcatList_prediction.source}, sum(selC), 1);
    C_eccodes    = repmat({[]}, sum(selC), 1);
    C_wc         = nan(sum(selC),1);
    C_orgn       = nan(sum(selC),1);

    % D) FUZZY (lower priority)
    D_rxns       = kcatList_fuzzy.rxns(selD);
    D_genes      = repmat({''}, sum(selD), 1);
    D_subs       = kcatList_fuzzy.substrates(selD);
    D_kcats      = kcatList_fuzzy.kcats(selD);
    D_source     = repmat({kcatList_fuzzy.source}, sum(selD), 1);
    D_eccodes    = getfield_or(kcatList_fuzzy, 'eccodes', numel(kcatList_fuzzy.rxns), '');
    D_eccodes    = D_eccodes(selD);
    D_wc         = wc(selD);
    D_orgn       = orgn(selD);

    % Concatenate in priority order A → B → C → D
    mergedKcatList.kcatSource  = [A_source;   B_source;   C_source;   D_source];
    mergedKcatList.rxns        = [A_rxns;     B_rxns;     C_rxns;     D_rxns];
    mergedKcatList.genes       = [A_genes;    B_genes;    C_genes;    D_genes];
    mergedKcatList.substrates  = [A_subs;     B_subs;     C_subs;     D_subs];
    mergedKcatList.kcats       = [A_kcats;    B_kcats;    C_kcats;    D_kcats];
    mergedKcatList.eccodes     = [A_eccodes;  B_eccodes;  C_eccodes;  D_eccodes];
    mergedKcatList.wildcardLvl = [A_wc;       B_wc;       C_wc;       D_wc];
    mergedKcatList.origin      = [A_orgn;     B_orgn;     C_orgn;     D_orgn];

    % Overall source tag
    mergedKcatList.source = 'Merged complete + fuzzy + prediction';
end
