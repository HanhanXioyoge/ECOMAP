function [model, geneIdx, enzIdx, protMetIdx] = addGeneEC(model, geneID, varargin)
% ADDGENEEC  Add a gene to an enzyme-constrained model, optionally adding its enzyme,
%            pseudo-metabolite (integrated), and a usage reaction (integrated, model-level only).
%
% Required:
%   model   : ecModel-like struct; must have model.enzymeConstraints.ecModeltype (non-empty)
%   geneID  : char/string, ID to append to model.genes (if not exists)
%
% Name-Value Options:
%   % ---- Gene layer ----
%   'GeneShortName'       : char/string, stored in model.geneShortNames (default: GeneName or geneID)
%   'GeneMiriam'          : struct or [], pushed into model.geneMiriams (default: [])
%   'OnDuplicate'         : 'skip'(default) | 'overwrite' | 'error'  % duplicate policy for GENE only
%
%   % ---- Enzyme layer ----
%   'AddEnzymeConstraint' : logical (default false). If true, also register in enzymeConstraints.*
%   'UniProtID'           : char/string, required when AddEnzymeConstraint=true (e.g., 'P0A799')
%   'MW'                  : double scalar; enzyme molecular weight (Da) (default NaN)
%   'Sequence'            : char/string; protein sequence (default '')
%   'PDB'                 : char/string; if empty, auto set to 'seq<index>.pdb' after insertion
%
%   % ---- Pseudo-metabolite / usage reaction (integrated only) ----
%   'ProtComp'            : compartment for 'prot_<UniProtID>' (index or ID like 'c'), default 'c'
%   'PoolMetID'           : char/string for protein pool metabolite (default 'prot_pool')
%   'PoolComp'            : compartment for PoolMetID (default: same as ProtComp)
%   'AddUsageRxn'         : logical (default true) — only effective in 'integrated'
%   'UsageOnDuplicate'    : 'skip'(default) | 'overwrite' | 'error'  % duplicate policy for usage rxn
%
% Returns:
%   model      : updated model
%   geneIdx    : index into model.genes
%   enzIdx     : index into model.enzymeConstraints.enzymes (if enzyme added or existed), else []
%   protMetIdx : index into model.mets for 'prot_<UniProtID>' (if integrated & enzyme added), else []
%
% Notes:
%   - NEVER edits enzymeConstraints.rxns here (usage reaction goes to model-level ONLY).
%   - Expands rxnGeneMat (#rxns x #genes) by adding zero-column when a NEW gene is appended.
%   - For 'integrated', requires addMetaboliteEC.m on path to create 'prot_<UniProtID>' / 'prot_pool'.

    % ----------------------------- Guards ----------------------------------
    assert(isfield(model, 'enzymeConstraints') && ...
           isfield(model.enzymeConstraints, 'ecModeltype') && ...
           ~isempty(model.enzymeConstraints.ecModeltype), ...
           'addGeneEC:NotEcModel: Missing/empty model.enzymeConstraints.ecModeltype.');

    if ~isfield(model,'genes') || ~iscell(model.genes), model.genes = cell(0,1); end
    if ~isfield(model,'rxns')  || ~iscell(model.rxns),  model.rxns  = cell(0,1); end
    if ~isfield(model,'S')     || isempty(model.S)
        model.S = sparse(numel(model.mets), numel(model.rxns));
    end

    % ------------------------- Parse inputs --------------------------------
    p = inputParser;
    p.FunctionName = 'addGeneEC';
    addRequired(p, 'model');
    addRequired(p, 'geneID',  @(x)ischar(x) || isstring(x));

    % Gene
    addParameter(p, 'GeneName',            '',        @(x)ischar(x) || isstring(x));
    addParameter(p, 'GeneShortName',       '',        @(x)ischar(x) || isstring(x));
    addParameter(p, 'GeneMiriam',          [],        @(x)isstruct(x) || isempty(x));
    addParameter(p, 'OnDuplicate',         'skip',    @(x)any(strcmpi(x,{'skip','overwrite','error'})));

    % Enzyme
    addParameter(p, 'AddEnzymeConstraint', false,     @(x)islogical(x) && isscalar(x));
    addParameter(p, 'UniProtID',           '',        @(x)ischar(x) || isstring(x));
    addParameter(p, 'MW',                  NaN,       @(x)isnumeric(x) && isscalar(x));
    addParameter(p, 'Sequence',            '',        @(x)ischar(x) || isstring(x));
    addParameter(p, 'PDB',                 '',        @(x)ischar(x) || isstring(x));

    % Integrated pseudo-metabolite & usage reaction
    addParameter(p, 'ProtComp',            'c',       @(x)ischar(x) || isstring(x) || (isnumeric(x)&&isscalar(x)));
    addParameter(p, 'PoolMetID',           'prot_pool', @(x)ischar(x) || isstring(x));
    addParameter(p, 'PoolComp',            '',        @(x)ischar(x) || isstring(x) || (isnumeric(x)&&isscalar(x)));
    addParameter(p, 'AddUsageRxn',         true,      @(x)islogical(x) && isscalar(x));
    addParameter(p, 'UsageOnDuplicate',    'skip',    @(x)any(strcmpi(x,{'skip','overwrite','error'})));

    parse(p, model, geneID, varargin{:});

    geneID     = char(p.Results.geneID);
    geneName   = char(p.Results.GeneName);
    geneSName  = char(p.Results.GeneShortName);
    geneMiriam = p.Results.GeneMiriam;
    dupGene    = lower(char(p.Results.OnDuplicate));

    addEnz     = p.Results.AddEnzymeConstraint;
    uniID      = char(p.Results.UniProtID);
    mw_in      = p.Results.MW;
    seq_in     = char(p.Results.Sequence);
    pdb_in     = char(p.Results.PDB);

    protComp   = p.Results.ProtComp;
    poolMetID  = char(p.Results.PoolMetID);
    poolComp   = p.Results.PoolComp;
    if isempty(poolComp), poolComp = protComp; end
    addUsage   = p.Results.AddUsageRxn;
    dupUsage   = lower(char(p.Results.UsageOnDuplicate));

    if addEnz && isempty(uniID)
        error('addGeneEC:UniProtMissing', 'UniProtID is required when AddEnzymeConstraint=true.');
    end

    % ----------------- Align gene-level arrays & matrices ------------------
    nGenesOld = numel(model.genes);
    nRxns     = numel(model.rxns);

    model = ensureField(model,'geneShortNames', @() repmat({''}, nGenesOld,1));
    model = ensureField(model,'geneMiriams',    @() cell(nGenesOld,1));

    if ~isfield(model,'rxnGeneMat') || isempty(model.rxnGeneMat)
        model.rxnGeneMat = sparse(nRxns, nGenesOld);
    else
        [rG,cG] = size(model.rxnGeneMat);
        if rG ~= nRxns || cG ~= nGenesOld
            error('addGeneEC:SizeMismatch', 'rxnGeneMat must be [#rxns x #genes].');
        end
    end

    % ---------------------- Add / reuse the gene --------------------------
    geneIdx = find(strcmp(model.genes, geneID), 1, 'first');
    if ~isempty(geneIdx)
        switch dupGene
            case 'skip'
                % leave untouched
            case 'overwrite'
                if ~isempty(geneSName), model.geneShortNames{geneIdx} = geneSName; end
                if ~isempty(geneMiriam),model.geneMiriams{geneIdx}    = geneMiriam;end
            case 'error'
                error('addGeneEC:DuplicateGene','Gene "%s" already exists.', geneID);
        end
    else
        if isempty(geneSName)
            geneSName = ternary(~isempty(geneName), geneName, geneID);
        end
        model.genes           = [model.genes;           {geneID}];
        model.geneShortNames  = [model.geneShortNames;  {geneSName}];
        model.geneMiriams     = [model.geneMiriams;     {geneMiriam}];
        model.rxnGeneMat      = [model.rxnGeneMat, sparse(nRxns,1)];  % add zero column
        geneIdx               = nGenesOld + 1;
    end

    % ---------------- EnzymeConstraints (enzymes/genes/mw/sequence/PDB) ----
    enzIdx     = [];
    protMetIdx = [];

    if addEnz
        EC = model.enzymeConstraints;
        if ~isfield(EC,'enzymes')   || isempty(EC.enzymes),   EC.enzymes   = cell(0,1); end
        if ~isfield(EC,'genes')     || isempty(EC.genes),     EC.genes     = cell(0,1); end
        if ~isfield(EC,'rxnEnzMat') || isempty(EC.rxnEnzMat), EC.rxnEnzMat = sparse(numel(EC.rxns), numel(EC.enzymes)); end

        nEnz = numel(EC.enzymes);
        EC = ensureECField(EC, 'mw',       nan(nEnz,1));      % double
        EC = ensureECField(EC, 'sequence', cell(nEnz,1));     % cellstr
        EC = ensureECField(EC, 'PDB',      cell(nEnz,1));     % cellstr

        if numel(EC.genes) ~= nEnz
            error('addGeneEC:ECPairMismatch', 'enzymeConstraints.genes and .enzymes must have same length.');
        end
        [rE,cE] = size(EC.rxnEnzMat);
        if rE ~= numel(EC.rxns) || cE ~= nEnz
            error('addGeneEC:ECRxnSz', 'enzymeConstraints.rxnEnzMat must be [#EC.rxns x #enzymes].');
        end

        % Duplicate ENZYME policy: ALWAYS SKIP & STOP downstream steps (as requested)
        hitE = find(strcmp(EC.enzymes, uniID), 1, 'first');
        if ~isempty(hitE)
            enzIdx = hitE;
            % do NOT modify EC.genes/mw/sequence/PDB; do NOT add prot_ metabolite; do NOT add usage rxn.
            model.enzymeConstraints = EC;
            return;
        end

        % Append new enzyme & paired gene (1:1 order); grow EC.rxnEnzMat by one column
        EC.enzymes   = [EC.enzymes;   {uniID}];
        EC.genes     = [EC.genes;     {geneID}];
        EC.rxnEnzMat = [EC.rxnEnzMat, sparse(size(EC.rxnEnzMat,1), 1)]; % add zero column

        nEnzNew      = nEnz + 1;
        EC = growECTo(EC, 'mw',       nEnzNew, nan);   EC.mw(nEnzNew)       = mw_in;
        EC = growECTo(EC, 'sequence', nEnzNew, {''});  EC.sequence{nEnzNew} = seq_in;
        EC = growECTo(EC, 'PDB',      nEnzNew, {''});
        if isempty(strtrim(pdb_in))
            EC.PDB{nEnzNew} = sprintf('seq%d.pdb', nEnzNew);
        else
            EC.PDB{nEnzNew} = pdb_in;
        end
        enzIdx = nEnzNew;

        % Commit EC update
        model.enzymeConstraints = EC;

        % ------------- integrated: add 'prot_<UniProtID>' and usage reaction -------------
        if strcmpi(model.enzymeConstraints.ecModeltype, 'integrated')
            % (1) ensure pseudo-metabolites exist
            assert(exist('addMetaboliteEC','file')==2 || exist('addMetaboliteEC','file')==6, ...
                   'addGeneEC:MissingAddMetaboliteEC: addMetaboliteEC.m must be on path for integrated models.');

            protMetID  = ['prot_', uniID];
            protMetIdx = find(strcmp(model.mets, protMetID), 1, 'first');
            if isempty(protMetIdx)
                model = addMetaboliteEC(model, protMetID, protComp, ...
                        'Name',   protMetID, ...
                        'Formula','', ...
                        'Charge', NaN, ...
                        'Notes',  'Enzyme-usage pseudometabolite');
                protMetIdx = numel(model.mets);
            end

            poolIdx = find(strcmp(model.mets, poolMetID), 1, 'first');
            if isempty(poolIdx)
                model = addMetaboliteEC(model, poolMetID, poolComp, ...
                        'Name',   poolMetID, ...
                        'Formula','', ...
                        'Charge', NaN, ...
                        'Notes',  'Total protein pool');
                poolIdx = numel(model.mets);
            end

            % (2) add usage reaction at MODEL level only (not EC layer)
            if addUsage
                usageID = ['usage_prot_', uniID];
                existRxn = find(strcmp(model.rxns, usageID), 1, 'first');
                switch dupUsage
                    case 'skip'
                        if ~isempty(existRxn)
                            % do nothing
                            return;
                        end
                        model = appendUsageRxn(model, usageID, protMetIdx, poolIdx);
                    case 'overwrite'
                        if ~isempty(existRxn)
                            model = overwriteUsageRxn(model, existRxn, protMetIdx, poolIdx);
                        else
                            model = appendUsageRxn(model, usageID, protMetIdx, poolIdx);
                        end
                    case 'error'
                        if ~isempty(existRxn)
                            error('addGeneEC:DuplicateUsageRxn','Usage reaction "%s" already exists.', usageID);
                        end
                        model = appendUsageRxn(model, usageID, protMetIdx, poolIdx);
                end
            end
        end
    end
end

% =========================== Local helpers =================================

function model = appendUsageRxn(model, usageID, protMetIdx, poolIdx)
    % Ensure rxn-level vectors aligned
    nRxns = numel(model.rxns);
    model = ensureField(model,'rxnNames',   @() repmat({''}, nRxns,1));
    model = ensureField(model,'lb',         @() zeros(nRxns,1));
    model = ensureField(model,'ub',         @() zeros(nRxns,1));
    model = ensureField(model,'rev',        @() false(nRxns,1));
    model = ensureField(model,'c',          @() zeros(nRxns,1));
    model = ensureField(model,'rxnMiriams', @() cell(nRxns,1));
    if ~isfield(model,'eccodes') || isempty(model.eccodes)
        model.eccodes = repmat({''}, nRxns,1);
    end
    % S: add new column
    nMets = numel(model.mets);
    if size(model.S,1) ~= nMets
        error('addGeneEC:SRows','S rows must equal numel(model.mets).');
    end
    model.S = [model.S, sparse(nMets,1)];
    newCol  = size(model.S,2);
    model.S(protMetIdx, newCol) = -1;
    model.S(poolIdx,    newCol) = +1;

    % Append rxn record
    model.rxns      = [model.rxns; {usageID}];
    model.rxnNames  = [model.rxnNames; {usageID}];
    model.lb        = [model.lb; -1000];
    model.ub        = [model.ub;  1];
    model.rev       = [model.rev; true];  % lb<0 & ub>0
    model.c         = [model.c;   0];
    model.rxnMiriams= [model.rxnMiriams; {struct('name',{},'value',{})}];
    model.eccodes   = [model.eccodes; {''}];

    % grow rxnGeneMat by one row (all zeros)
    if ~isfield(model,'rxnGeneMat') || isempty(model.rxnGeneMat)
        model.rxnGeneMat = sparse(numel(model.rxns)-1, numel(model.genes));
    end
    [rG,cG] = size(model.rxnGeneMat);
    if rG ~= numel(model.rxns)-1
        if rG < numel(model.rxns)-1
            model.rxnGeneMat(numel(model.rxns)-1,1) = 0;
        else
            error('addGeneEC:rxnGeneMatRows','rxnGeneMat rows inconsistent.');
        end
    end
    if cG < numel(model.genes)
        model.rxnGeneMat(:, cG+1:numel(model.genes)) = 0;
    end
    model.rxnGeneMat = [model.rxnGeneMat; sparse(1, numel(model.genes))];

    % ensure grRules length
    model = ensureField(model,'grRules', @() repmat({''}, numel(model.rxns)-1,1));
    model.grRules = [model.grRules; {''}];
end

function model = overwriteUsageRxn(model, rxnIdx, protMetIdx, poolIdx)
    % Overwrite stoich & bounds & names in-place
    model.S(:, rxnIdx) = 0;
    model.S(protMetIdx, rxnIdx) = -1;
    model.S(poolIdx,    rxnIdx) = +1;

    model.rxnNames{rxnIdx} = model.rxns{rxnIdx};
    model.lb(rxnIdx)  = -1000;
    model.ub(rxnIdx)  = 1;
    model.rev(rxnIdx) = true;
    model.c(rxnIdx)   = 0;
    model.rxnMiriams{rxnIdx} = struct('name',{},'value',{});
    if isfield(model,'eccodes') && numel(model.eccodes)>=rxnIdx
        model.eccodes{rxnIdx} = '';
    end
    % keep rxnGeneMat row as-is (all zeros typically)
    % keep grRules as-is (empty)
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function model = ensureField(model, fld, initFcn)
    if ~isfield(model, fld) || isempty(model.(fld))
        model.(fld) = initFcn();
    end
end

function EC = ensureECField(EC, field, initVal)
    n = numel(EC.enzymes);
    if ~isfield(EC, field) || isempty(EC.(field))
        EC.(field) = initVal;
        if iscell(EC.(field))
            EC.(field) = ensureCellLen(EC.(field), n, '');
        else
            EC.(field) = ensureNumLen(EC.(field),  n, NaN);
        end
    else
        if iscell(EC.(field))
            if numel(EC.(field)) ~= n, error('addGeneEC:LenMismatchECField','enzymeConstraints.%s length mismatch.',field); end
        end
    end
end

function EC = growECTo(EC, field, N, fillVal)
    if iscell(EC.(field))
        cur = numel(EC.(field));
        if cur < N
            EC.(field){N} = [];
            for k = cur+1:N, EC.(field){k} = fillVal{:}; end
        end
    else
        cur = numel(EC.(field));
        if cur < N
            EC.(field)(N,1) = fillVal;
            if cur+1 <= N, EC.(field)(cur+1:N,1) = fillVal; end
        end
    end
end

function c = ensureCellLen(c, n, fillVal)
    if isempty(c), c = cell(0,1); end
    cur = numel(c);
    if cur < n
        c{n} = [];
        for k = cur+1:n, c{k} = fillVal; end
    elseif cur > n
        c = c(1:n);
    end
end

function x = ensureNumLen(x, n, fillVal)
    if isempty(x), x = []; end
    cur = numel(x);
    if cur < n
        x(n,1) = fillVal;
        if cur+1 <= n, x(cur+1:n,1) = fillVal; end
    elseif cur > n
        x = x(1:n);
    end
end
