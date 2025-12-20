function model = buildBiomassPseudoRxns(model, comps, compTypes, bioRxnID)
% buildBiomassPseudoRxns
% -------------------------------------------------------------------------
% Build biomass composition pseudo-reactions (protein / carbohydrate / RNA /
% DNA / lipid) from an existing biomass reaction, using a component list.
%
% For each requested compType in compTypes:
%   - Collect the biomass stoichiometry of components of that type
%     (from the biomass reaction bioRxnID).
%   - Create a pseudo-metabolite (e.g. "protein") as a biomass pool.
%   - Create a pseudo-reaction (e.g. "protein pseudoreaction"):
%
%         sum_i (coeff_i * component_i)  -->  protein_pool
%
%   - The stoichiometric coefficients for component_i are copied from the
%     biomass reaction, so that GECKO-like functions (sumBioMass, etc.)
%     can compute mass fractions from these pseudo-reactions.
%
% INPUT:
%   model     : ecModel-like structure (must be compatible with
%               addMetaboliteEC / addReactionEC, i.e. have
%               model.enzymeConstraints.ecModeltype).
%   comps     : N x 3 cell array:
%                   { metID, MW, classLetter }
%               Example (E. coli amino acids, protein "P" class):
%                   { 'ala__L_c',  89.09, 'P';
%                     'cys__L_c', 121.16, 'P';
%                     ... }
%               classLetter ∈ {'P','C','R','D','L'}.
%   compTypes : cellstr or char, subset of {'P','C','R','D','L'}
%               e.g. {'P'} or {'P','C','R','D','L'}
%   bioRxnID  : biomass reaction ID in model.rxns, e.g.
%               'BIOMASS_Ec_iML1515_core_75p37M'
%
% OUTPUT:
%   model     : updated model with additional pseudo-metabolites and
%               pseudo-reactions for each requested compType.
%
% NOTES:
%   - This function does NOT modify the biomass reaction itself; it only
%     reads its stoichiometry to construct composition pseudo-reactions.
%   - Pseudo-reaction naming scheme:
%         compType='P' -> 'protein pseudoreaction'
%         compType='C' -> 'carbohydrate pseudoreaction'
%         compType='R' -> 'RNA pseudoreaction'
%         compType='D' -> 'DNA pseudoreaction'
%         compType='L' -> 'lipid backbone pseudoreaction'
%     Pseudo-metabolite names: 'protein', 'carbohydrate', 'RNA', 'DNA',
%     'lipid backbone'; IDs like 'protein_biomass_c'.
%   - Pseudo-metabolites are placed in the same compartment as their
%     components (inferred from metComps), or in 'c' if unavailable.
% -------------------------------------------------------------------------

    % -------- Basic checks --------
    if nargin < 3 || isempty(compTypes)
        error('compTypes (e.g. {''P''} or {''P'',''C''}) must be provided.');
    end
    if nargin < 4 || isempty(bioRxnID)
        error('bioRxnID (biomass reaction ID) must be provided.');
    end

    if ischar(compTypes) || isstring(compTypes)
        compTypes = {char(compTypes)};
    end

    if ~iscell(comps) || size(comps,2) < 3
        error('comps must be an N x 3 cell array: {metID, MW, classLetter}.');
    end

    if ~isfield(model,'S') || isempty(model.S)
        error('buildBiomassPseudoRxns:NoS', ...
              'Model.S (stoichiometric matrix) is missing or empty.');
    end
    if ~isfield(model,'rxns')
        error('buildBiomassPseudoRxns:NoRxns','model.rxns must exist.');
    end
    if ~isfield(model,'mets')
        error('buildBiomassPseudoRxns:NoMets','model.mets must exist.');
    end

    % Biomass reaction index
    bioIdx = find(strcmp(model.rxns, bioRxnID));
    if isempty(bioIdx)
        error('Biomass reaction "%s" not found in model.rxns.', bioRxnID);
    elseif numel(bioIdx) > 1
        error('Multiple reactions match biomass ID "%s".', bioRxnID);
    end

    % For convenience
    Sbio  = model.S(:, bioIdx);   % stoichiometric column of biomass reaction

    % -------- Loop over component types (P/C/R/D/L) --------
    for k = 1:numel(compTypes)
        compType = upper(compTypes{k});  % 'P','C','R','D','L'

        % Find rows in 'comps' for this compType
        classCol = comps(:,3);
        maskType = cellfun(@(x)strcmpi(x,compType), classCol);
        if ~any(maskType)
            warning('No components with class "%s" found in comps. Skipping.', compType);
            continue;
        end

        compMetIDs = comps(maskType, 1);  % metIDs (e.g. 'ala__L_c', etc.)

        % Map these component metIDs to indices in model.mets
        compMetIdxs = zeros(numel(compMetIDs), 1);
        for i = 1:numel(compMetIDs)
            metID = char(compMetIDs{i});
            idx = find(strcmp(model.mets, metID), 1, 'first');
            if isempty(idx)
                error('Component metabolite "%s" (class "%s") not found in model.mets.', ...
                      metID, compType);
            end
            compMetIdxs(i) = idx;
        end

        % Determine compartment for the pseudo-metabolite:
        %   - If metComps exists, use the compartment of the first component.
        %   - Otherwise, default to 'c'.
        if isfield(model,'metComps') && ~isempty(model.metComps)
            baseCompIdx = model.metComps(compMetIdxs(1));
        else
            baseCompIdx = 'c';  % will be resolved by addMetaboliteEC
        end

        % Get pseudo-metabolite & pseudo-reaction names/IDs
        [pseudoMetName, pseudoMetID, pseudoRxnName, pseudoRxnID] = ...
            getCompPseudoNames(compType);

        % Ensure pseudo-metabolite exists (search by metID first, then by metName)
        pseudoMetIdx = [];
        hitID  = find(strcmp(model.mets, pseudoMetID), 1, 'first');
        if ~isempty(hitID)
            pseudoMetIdx = hitID;
        elseif isfield(model,'metNames')
            hitName = find(strcmp(model.metNames, pseudoMetName), 1, 'first');
            if ~isempty(hitName)
                pseudoMetIdx = hitName;
            end
        end

        if isempty(pseudoMetIdx)
            % Create pseudo-metabolite via addMetaboliteEC
            [model, pseudoMetIdx] = addMetaboliteEC(model, pseudoMetID, baseCompIdx, ...
                                        'Name', pseudoMetName);
            % After adding, model.mets and model.S are extended by one row
            Sbio  = model.S(:, bioIdx);   % refresh Sbio reference
        end

        % If pseudo-reaction already exists, skip (so function is idempotent)
        if isfield(model,'rxnNames') && any(strcmp(model.rxnNames, pseudoRxnName))
            warning('Pseudo-reaction "%s" already exists. Skipping.', pseudoRxnName);
            continue;
        end
        if any(strcmp(model.rxns, pseudoRxnID))
            warning('Pseudo-reaction ID "%s" already exists. Skipping.', pseudoRxnID);
            continue;
        end

        % -------- Build stoichiometry for the new pseudo-reaction --------
        % We copy the biomass reaction stoichiometry for the relevant
        % component metabolites, and add the pseudo-metabolite as product.
        metsRxn   = cell(numel(compMetIdxs)+1, 1);
        stoichRxn = zeros(numel(compMetIdxs)+1, 1);

        for j = 1:numel(compMetIdxs)
            idx = compMetIdxs(j);
            coeff = Sbio(idx);   % stoichiometric coefficient in biomass reaction
            if coeff == 0
                % This component metabolite does not appear in biomass reaction,
                % which is unexpected if comps is consistent with biomass.
                warning('Metabolite "%s" (class "%s") has zero coeff. in biomass; skipping it.', ...
                        model.mets{idx}, compType);
                continue;
            end
            metsRxn{j}   = model.mets{idx};
            stoichRxn(j) = coeff;
        end

        % Product: pseudo-metabolite (e.g., protein_biomass_c) with coeff +1
        metsRxn{end}   = model.mets{pseudoMetIdx};  % can also use pseudoMetID
        stoichRxn(end) = 1;

        % Remove any empty entries that might come from skipped components
        emptyMask = cellfun(@isempty, metsRxn);
        metsRxn(emptyMask)   = [];
        stoichRxn(emptyMask) = [];

        if isempty(metsRxn)
            warning('No non-zero biomass components for class "%s". Skipping pseudo-reaction.', compType);
            continue;
        end

        % -------- Add pseudo-reaction via addReactionEC --------
        % Use irreversible forward direction; bounds can be adjusted as needed.
        lb = 0;
        ub = 1000;

        % We do not specify genes or enzymes here; this is a purely
        % compositional pseudo-reaction.
        model = addReactionEC(model, pseudoRxnID, metsRxn, stoichRxn, lb, ub, ...
                              'RxnName',   pseudoRxnName, ...
                              'OnDuplicate','skip');

        model.S(compMetIdxs, bioIdx) = 0;
        model.S(pseudoMetIdx, bioIdx) = -1;
        nMets = size(model.S, 1);     % Number of metabolites = number of rows of S
        model.csense = repmat('E', nMets, 1);
    end
end

function [metName, metID, rxnName, rxnID] = getCompPseudoNames(compType)
% getCompPseudoNames
% -------------------------------------------------------------------------
% Map compType ('P','C','R','D','L') to pseudo-metabolite and
% pseudo-reaction names/IDs.
%
% OUTPUT:
%   metName : pseudo-metabolite full name (e.g. 'protein')
%   metID   : pseudo-metabolite ID (e.g. 'protein_biomass_c')
%   rxnName : pseudo-reaction name  (e.g. 'protein pseudoreaction')
%   rxnID   : pseudo-reaction ID    (e.g. 'protein_pseudoreaction')
% -------------------------------------------------------------------------

    switch upper(compType)
        case 'P'
            metName = 'protein';
            metID   = 'protein_biomass_c';
        case 'C'
            metName = 'carbohydrate';
            metID   = 'carbohydrate_biomass_c';
        case 'R'
            metName = 'RNA';
            metID   = 'RNA_biomass_c';
        case 'D'
            metName = 'DNA';
            metID   = 'DNA_biomass_c';
        case 'L'
            metName = 'lipid backbone';
            metID   = 'lipid_biomass_c';
        otherwise
            error('Unknown compType "%s". Must be one of P/C/R/D/L.', compType);
    end

    rxnName = [metName ' pseudoreaction'];          % e.g. 'protein pseudoreaction'
    rxnID   = lower(regexprep(rxnName,'\s+','_'));  % e.g. 'protein_pseudoreaction'
end