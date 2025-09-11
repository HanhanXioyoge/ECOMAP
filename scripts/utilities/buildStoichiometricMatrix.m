function S = buildStoichiometricMatrix(mets, reactions)
% BUILDSTOICHIOMETRICMATRIX Construct a stoichiometric matrix 
% from the structure of the reactive metabolite
% input：
%   mets       - Metabolite ID list (array of cells)
%   reactions  - An array of reaction structures, each element 
%                contains a metabolite substructure
% output：
%   S          - sparse stoichiometric matrix（mets × reactions）

    numMets = numel(mets);
    numRxns = numel(reactions);
    
    rowIdx = [];
    colIdx = [];
    values = [];
    
    for rxnIdx = 1:numRxns
        metabolites = reactions{rxnIdx}.metabolites;
        
        metIDs = fieldnames(metabolites);
        stoich = struct2array(metabolites); 
        
        if ~isnumeric(stoich)
            error('The stoichiometric number of %d of the reaction should be of numerical type', rxnIdx);
        end
        
        [~, metPositions] = ismember(metIDs, mets);
        
        valid = metPositions ~= 0;
        metPositions = metPositions(valid);
        stoich = stoich(valid);
        
        rowIdx = [rowIdx; metPositions];
        colIdx = [colIdx; rxnIdx * ones(length(metPositions), 1)];
        values = [values; stoich(:)]; 
    end
    
    S = sparse(rowIdx, colIdx, values, numMets, numRxns);
    
    % dimensionality verification
    assert(size(S,1) == numMets, 'Inconsistent metabolite dimensions: expected %d actual %d', numMets, size(S,1));
    assert(size(S,2) == numRxns, 'Inconsistent response dimensions: expected %d actual %d', numRxns, size(S,2));
end