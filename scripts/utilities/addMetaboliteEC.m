function [model, newMetIdx] = addMetaboliteEC(model, metID, comp, varargin)
% ADDMETABOLITEC  Add a metabolite row into an enzyme-constrained model.
%
% Usage:
%   model = addMetaboliteEC(model, metID, comp)
%   model = addMetaboliteEC(model, metID, comp, 'Name', 'ATP', ...
%                           'Formula','C10H16N5O13P3', 'Charge', -4, 'Notes', 'new metabolite')
%   [model, idx] = addMetaboliteEC(...)
%
% Required:
%   model  : ecModel-like struct
%   metID  : char/string, metabolite abbreviation to append to model.mets
%   comp   : (a) numeric index into model.comps
%            (b) char/string matching model.comps or model.compNames
%
% Name-Value Pairs (optional):
%   'Name'       : metabolite full name   (default: '')
%   'Formula'    : chemical formula       (default: '')
%   'Charge'     : numeric charge         (default: NaN)
%   'Notes'      : free text notes        (default: '')
%   'Smiles'     : metabolite Smiles      (default: '')
%   'InChIKey'   : metabolite InChIKey    (default: '')
%   'MetaNetXID' : metabolite MetaNetXID  (default: '')
%
% Behavior:
%   - Verifies model.enzymeConstraints.ecModeltype is present and nonempty.
%   - Ensures all targeted met-* vectors are length-aligned after append:
%       mets, metNames, metComps, metFormulas, metCharges (or metCharge),
%       metNotes, metMiriams, b, and (if present) metSmiles / metInChIKey /
%       metMetaNetXID.
%   - If model does NOT have metSmiles/metInChIKey/metMetaNetXID:
%       * these fields are NOT created;
%       * giving a non-empty Smiles/InChIKey/MetaNetXID will raise an error.
%   - If model DOES have those fields:
%       * they must be aligned (length == numel(mets));
%       * if no value is provided, a default '' is appended.
%   - Expands S by adding a zero row (sparse). If S missing, it will be created.
%   - Does NOT connect the new metabolite to any reaction (all zeros in S).
%
% Returns:
%   model     : updated model
%   newMetIdx : index of the newly appended metabolite (optional)

    % ----------------------------- Guard rails -----------------------------
    if ~isfield(model, 'enzymeConstraints') || ...
       ~isfield(model.enzymeConstraints, 'ecModeltype') || ...
       isempty(model.enzymeConstraints.ecModeltype)
        error('addMetaboliteEC:NotEcModel', ...
            ['This model does not appear to be enzyme-constrained: ', ...
             'missing or empty model.enzymeConstraints.ecModeltype.']);
    end

    if ~isfield(model,'mets') || ~iscell(model.mets)
        error('addMetaboliteEC:BadModel', 'model.mets must exist and be a cell array.');
    end
    if ~isfield(model,'metComps')
        model.metComps = zeros(numel(model.mets),1);
    end

    % ------------------------- Parse inputs & defaults ---------------------
    p = inputParser;
    p.FunctionName = 'addMetaboliteEC';
    addRequired(p, 'model');
    addRequired(p, 'metID',  @(x)ischar(x) || isstring(x));
    addRequired(p, 'comp');
    addParameter(p, 'Name',       '', @(x)ischar(x) || isstring(x));
    addParameter(p, 'Formula',    '', @(x)ischar(x) || isstring(x));
    addParameter(p, 'Charge',    NaN, @(x)isnumeric(x) && isscalar(x));
    addParameter(p, 'Notes',      '', @(x)ischar(x) || isstring(x));
    addParameter(p, 'Smiles',     '', @(x)ischar(x) || isstring(x));
    addParameter(p, 'InChIKey',   '', @(x)ischar(x) || isstring(x));
    addParameter(p, 'MetaNetXID', '', @(x)ischar(x) || isstring(x));
    parse(p, model, metID, comp, varargin{:});

    metID       = char(p.Results.metID);
    metName     = char(p.Results.Name);
    metForm     = char(p.Results.Formula);
    metChg      = p.Results.Charge;
    metNote     = char(p.Results.Notes);
    Smiles      = char(p.Results.Smiles);
    InChIKey    = char(p.Results.InChIKey);
    MetaNetXID  = char(p.Results.MetaNetXID);

    % -------------------- Resolve compartment to numeric index -------------
    compIdx = resolveCompIndex(model, comp);

    % ------------------ Prepare target met-* fields to align ---------------
    nOld = numel(model.mets);

    if ~isfield(model,'metNames') || isempty(model.metNames)
        model.metNames = repmat({''}, nOld, 1);
    elseif numel(model.metNames) ~= nOld
        error('addMetaboliteEC:LengthMismatch', 'model.metNames length must match model.mets.');
    end

    if ~isfield(model,'metFormulas') || isempty(model.metFormulas)
        model.metFormulas = repmat({''}, nOld, 1);
    elseif numel(model.metFormulas) ~= nOld
        error('addMetaboliteEC:LengthMismatch', 'model.metFormulas length must match model.mets.');
    end

    if ~isfield(model,'metNotes') || isempty(model.metNotes)
        model.metNotes = repmat({''}, nOld, 1);
    elseif numel(model.metNotes) ~= nOld
        error('addMetaboliteEC:LengthMismatch', 'model.metNotes length must match model.mets.');
    end

    % metCharges / metCharge compatibility
    hasPlural   = isfield(model,'metCharges') && ~isempty(model.metCharges);
    hasSingular = isfield(model,'metCharge')  && ~isempty(model.metCharge);
    if ~(hasPlural || hasSingular)
        model.metCharges = nan(nOld,1);
        hasPlural = true;
    end
    if hasPlural   && numel(model.metCharges) ~= nOld
        error('addMetaboliteEC:LengthMismatch','model.metCharges length must match model.mets.');
    end
    if hasSingular && numel(model.metCharge)  ~= nOld
        error('addMetaboliteEC:LengthMismatch','model.metCharge length must match model.mets.');
    end

    if numel(model.metComps) ~= nOld
        error('addMetaboliteEC:LengthMismatch', 'model.metComps length must match model.mets.');
    end
    
    % ---------------- Smiles / InChIKey / MetaNetXID  -------------
    hasSmilesField     = isfield(model,'metSmiles');
    hasInChIKeyField   = isfield(model,'metInChIKey');
    hasMetaNetXIDField = isfield(model,'metMetaNetXID');

    % metSmiles
    if hasSmilesField
        if ~isempty(model.metSmiles) && numel(model.metSmiles) ~= nOld
            error('addMetaboliteEC:LengthMismatch', 'model.metSmiles length must match model.mets.');
        end
    else
        if ~isempty(strtrim(Smiles))
            error('addMetaboliteEC:FieldMissing', ...
                  'Model has no metSmiles field but a Smiles value was provided.');
        end
    end

    % metInChIKey
    if hasInChIKeyField
        if ~isempty(model.metInChIKey) && numel(model.metInChIKey) ~= nOld
            error('addMetaboliteEC:LengthMismatch', 'model.metInChIKey length must match model.mets.');
        end
    else
        if ~isempty(strtrim(InChIKey))
            error('addMetaboliteEC:FieldMissing', ...
                  'Model has no metInChIKey field but an InChIKey value was provided.');
        end
    end

    % metMetaNetXID
    if hasMetaNetXIDField
        if ~isempty(model.metMetaNetXID) && numel(model.metMetaNetXID) ~= nOld
            error('addMetaboliteEC:LengthMismatch', 'model.metMetaNetXID length must match model.mets.');
        end
    else
        if ~isempty(strtrim(MetaNetXID))
            error('addMetaboliteEC:FieldMissing', ...
                  'Model has no metMetaNetXID field but a MetaNetXID value was provided.');
        end
    end

    % ---------------- NEW: Ensure metMiriams aligned -----------------------
    if ~isfield(model,'metMiriams') || isempty(model.metMiriams)
        model.metMiriams = cell(nOld,1);
    elseif numel(model.metMiriams) ~= nOld
        error('addMetaboliteEC:LengthMismatch', 'model.metMiriams length must match model.mets.');
    end

    if startsWith(string(metID), "prot_")
        emptyMiriam = struct('name', {{'sbo'}}, 'value', {{'SBO:0000252'}});
    else
        emptyMiriam = struct('name',{},'value',{});
    end

    % ---------------- NEW: Ensure b exists & aligned (m x 1) ---------------
    if ~isfield(model,'b') || isempty(model.b)
        model.b = zeros(nOld,1);
    elseif numel(model.b) ~= nOld
        error('addMetaboliteEC:LengthMismatch', 'model.b length must match model.mets (rows of S).');
    end

    % --------------------------- Append new entries ------------------------
    model.mets        = [model.mets;        {metID}];
    model.metNames    = [model.metNames;    {metName}];
    model.metFormulas = [model.metFormulas; {metForm}];
    model.metNotes    = [model.metNotes;    {metNote}];
    model.metComps    = [model.metComps;     compIdx];

    if hasPlural,   model.metCharges = [model.metCharges; metChg]; end
    if hasSingular, model.metCharge  = [model.metCharge;  metChg]; end

    % --------- Smiles / InChIKey / MetaNetXID -------------------
    if hasSmilesField
        if isempty(model.metSmiles)
            model.metSmiles = repmat({''}, nOld, 1);
        end
        if isempty(Smiles)
            model.metSmiles = [model.metSmiles; {''}];
        else
            model.metSmiles = [model.metSmiles; {Smiles}];
        end
    end

    if hasInChIKeyField
        if isempty(model.metInChIKey)
            model.metInChIKey = repmat({''}, nOld, 1);
        end
        if isempty(InChIKey)
            model.metInChIKey = [model.metInChIKey; {''}];
        else
            model.metInChIKey = [model.metInChIKey; {InChIKey}];
        end
    end

    if hasMetaNetXIDField
        if isempty(model.metMetaNetXID)
            model.metMetaNetXID = repmat({''}, nOld, 1);
        end
        if isempty(MetaNetXID)
            model.metMetaNetXID = [model.metMetaNetXID; {''}];
        else
            model.metMetaNetXID = [model.metMetaNetXID; {MetaNetXID}];
        end
    end

    % ---------------- NEW: Append metMiriams & b ---------------------------
    model.metMiriams  = [model.metMiriams; {emptyMiriam}];
    model.b           = [model.b; 0];

    newMetIdx = nOld + 1;

    % ---------------------------- Expand S matrix --------------------------
    if isfield(model,'S') && ~isempty(model.S)
        [r, c] = size(model.S);
        if r ~= nOld
            error('addMetaboliteEC:SRowsMismatch', ...
                'Size mismatch: size(model.S,1) = %d but numel(model.mets)-1 = %d.', r, nOld);
        end
        model.S = [model.S; sparse(1, c)];
    else
        nRxns = (isfield(model,'rxns') && ~isempty(model.rxns)) * numel(model.rxns);
        model.S = sparse(nOld+1, nRxns);
    end

    % --------------------- (Optional) Sanity post-checks -------------------
    nNow = numel(model.mets);
    assert(numel(model.metNames)    == nNow);
    assert(numel(model.metFormulas) == nNow);
    assert(numel(model.metNotes)    == nNow);
    assert(numel(model.metComps)    == nNow);
    if hasPlural,   assert(numel(model.metCharges) == nNow); end
    if hasSingular, assert(numel(model.metCharge)  == nNow); end
    assert(numel(model.metMiriams) == nNow);
    assert(numel(model.b)          == nNow);
    if hasSmilesField     && ~isempty(model.metSmiles)
        assert(numel(model.metSmiles) == nNow);
    end
    if hasInChIKeyField   && ~isempty(model.metInChIKey)
        assert(numel(model.metInChIKey) == nNow);
    end
    if hasMetaNetXIDField && ~isempty(model.metMetaNetXID)
        assert(numel(model.metMetaNetXID) == nNow);
    end
end

% ----------------------------- Helper -------------------------------------
function compIdx = resolveCompIndex(model, comp)
    if isnumeric(comp)
        compIdx = comp;
        if ~isscalar(compIdx) || compIdx < 1 || ...
           (isfield(model,'comps') && compIdx > numel(model.comps))
            error('addMetaboliteEC:BadCompIndex', 'Invalid compartment index.');
        end
        compIdx = double(compIdx);
        return;
    end
    if ~(ischar(comp) || isstring(comp))
        error('addMetaboliteEC:BadCompInput', 'comp must be a numeric index or a string.');
    end
    compStr = char(comp);
    if isfield(model,'comps') && ~isempty(model.comps)
        hit = find(strcmp(model.comps, compStr), 1, 'first');
        if ~isempty(hit), compIdx = hit; return; end
    end
    if isfield(model,'compNames') && ~isempty(model.compNames)
        hit = find(strcmp(model.compNames, compStr), 1, 'first');
        if ~isempty(hit), compIdx = hit; return; end
    end
    error('addMetaboliteEC:CompNotFound', ...
          'Compartment ''%s'' not found in model.comps or model.compNames.', compStr);
end
