function ecModel = cleanupUnusedEnzymes(ecModel)
% CLEANUPUNUSEDENZYMES
%   Remove enzymes (columns in rxnEnzMat) that are not used by any reaction
%   and clean up the corresponding gene / prot_* usage in the model.
%
%   Logic:
%   1) In ecModel.enzymeConstraints.rxnEnzMat, find enzyme columns that are
%      all zero (no reaction uses this enzyme).
%   2) Remove these columns from rxnEnzMat and from all enzyme-centric
%      fields (genes, enzymes, mw, sequence, PDB, geneIdx, enz2protIdx).
%   3) For the corresponding genes, if they are not used in any reaction
%      in model.rxnGeneMat, remove them from the COBRA model using
%      removemodelGenes.
%   4) For integrated ecModels that have protMets / enz2protIdx:
%      - Identify prot_* metabolites that belong ONLY to the removed
%        enzymes and are no longer referenced by any remaining enzyme.
%      - Remove their usage reactions (usage_prot_*) and, if prot_* is
%        not participating in any other reaction, remove the prot_* mets.
%      - Update protMets and enz2protIdx index mapping accordingly.

    % --- Basic guards ---
    if ~isfield(ecModel, 'enzymeConstraints') || ...
       ~isfield(ecModel.enzymeConstraints, 'rxnEnzMat')
        error('cleanupUnusedEnzymes:NoEC', ...
            'ecModel.enzymeConstraints.rxnEnzMat is missing.');
    end

    EC = ecModel.enzymeConstraints;

    % --- Check rxnEnzMat ---
    rxnEnzMat = EC.rxnEnzMat;
    if isempty(rxnEnzMat)
        % Nothing to do if there is no enzyme-reaction incidence matrix
        return;
    end

    % ---------------------------------------------------------------------
    % 1. Find enzyme columns that are never used (all-zero columns)
    % ---------------------------------------------------------------------
    % unusedCols: indices of enzyme columns with no non-zero entries.
    unusedCols = find(~any(rxnEnzMat ~= 0, 1));
    if isempty(unusedCols)
        % Every enzyme column is used by at least one reaction: nothing to clean
        return;
    end

    % Record genes & enzymes for those unused columns (for model-level cleanup)
    if isfield(EC, 'genes')
        genesToCheck = EC.genes(unusedCols);
    else
        genesToCheck = {};
    end
    if isfield(EC, 'enzymes')
        enzID_toCheck = EC.enzymes(unusedCols);
    else
        enzID_toCheck = {};
    end

    protIdxRemoved = [];

    % Build logical mask of columns to keep
    keepCols = true(1, size(rxnEnzMat, 2));
    keepCols(unusedCols) = false;

    % ---------------------------------------------------------------------
    % 2. Remove unused enzyme columns from enzymeConstraints
    % ---------------------------------------------------------------------
    % 2.1 rxnEnzMat: drop unused columns
    EC.rxnEnzMat = EC.rxnEnzMat(:, keepCols);

    % 2.2 All enzyme-centric arrays must be sliced along the same columns
    %     (genes / enzymes / mw / sequence / PDB / geneIdx / enz2protIdx)
    geneFields = {'genes','enzymes','mw','sequence','PDB','geneIdx','enz2protIdx'};
    for f = 1:numel(geneFields)
        fld = geneFields{f};
        if isfield(EC, fld) && ~isempty(EC.(fld))
            EC.(fld) = EC.(fld)(keepCols);
        end
    end

    % Write back to model
    ecModel.enzymeConstraints = EC;
    EC = ecModel.enzymeConstraints; % refresh local copy

    % ---------------------------------------------------------------------
    % 3. Remove genes from the COBRA model if they are not used anymore
    % ---------------------------------------------------------------------
    genesToRemove = {};
    for i = 1:numel(genesToCheck)
        gName = genesToCheck{i};
        if isempty(gName)
            continue;
        end

        % Find this gene in the COBRA model
        idxGene = find(strcmp(ecModel.genes, gName), 1);
        if isempty(idxGene)
            % Gene might already be removed; skip
            continue;
        end

        % Check if gene is still used in any reaction
        if isfield(ecModel, 'rxnGeneMat') && ~isempty(ecModel.rxnGeneMat)
            if ~any(ecModel.rxnGeneMat(:, idxGene))
                % Completely unused gene -> mark for removal
                genesToRemove{end+1} = gName; %#ok<AGROW>
            end
        end
    end

    % Remove genes from the model using RAVEN/COBRA helper if available
    if ~isempty(genesToRemove)
        if exist('removemodelGenes', 'file') == 2
            ecModel = removemodelGenes(ecModel, genesToRemove);
        else
            warning('cleanupUnusedEnzymes:NoRemovemodelGenes', ...
                ['removemodelGenes.m not found; genes in model were NOT removed. ', ...
                 'You may want to remove them manually: %s'], strjoin(genesToRemove, ', '));
        end
    end

    % ---------------------------------------------------------------------
    % 4. Remove usage reactions and prot_* metabolites for removed enzymes
    % ---------------------------------------------------------------------
    % We operate purely based on the UniProt IDs:
    %   prot ID   = ['prot_' UniProtID]
    %   usage ID  = ['usage_' prot ID]  -> e.g. 'usage_prot_P0A7Z4'
    %
    % For each removed UniProt ID:
    %   - If this UniProt ID is still present in EC.enzymes AFTER the column
    %     removal, we keep its usage reaction/metabolite (still used).
    %   - Only if no remaining enzyme uses this UniProt ID, we delete its
    %     usage reaction and prot_* metabolite (if disconnected).

    if isfield(EC, 'enzymes') && ~isempty(enzID_toCheck)
        % Remove empty entries and get unique UniProt IDs among removed enzymes
        removedEnzIDs = enzID_toCheck(~cellfun(@isempty, enzID_toCheck));
        if ~isempty(removedEnzIDs)
            uniqueRemoved = unique(removedEnzIDs);

            for k = 1:numel(uniqueRemoved)
                uniID = uniqueRemoved{k};

                % If this UniProt ID is still present in remaining EC.enzymes,
                % then some enzyme still uses it -> do NOT delete its usage.
                if any(strcmp(EC.enzymes, uniID))
                    continue;
                end

                % Build prot_* metabolite ID and usage_* reaction ID
                protID  = ['prot_'  uniID];       % e.g. 'prot_P0A7Z4'
                usageID = ['usage_' protID];      % e.g. 'usage_prot_P0A7Z4'

                % Remove usage reaction if present
                if ismember(usageID, ecModel.rxns)
                    try
                        ecModel = removeRxns(ecModel, usageID);
                    catch
                        warning('cleanupUnusedEnzymes:removeRxns', ...
                            'Failed to remove usage reaction %s. Please check removeRxns signature.', usageID);
                    end
                end

                % Remove prot_* metabolite if it is no longer connected
                metIdx = find(strcmp(ecModel.mets, protID), 1);
                if ~isempty(metIdx)
                    if all(ecModel.S(metIdx, :) == 0)
                        try
                            ecModel = removeMets(ecModel, protID);
                        catch
                            warning('cleanupUnusedEnzymes:removeMets', ...
                                'Failed to remove metabolite %s. Please check removeMets signature.', protID);
                        end
                    end
                end
            end
        end
    end
end
