function rxnGeneMat = buildRxnGeneMat(genes, grRules)
% BUILDRXNGENEMAT Constructs the reaction-gene association matrix from gene rules.
% Inputs:
%   genes   - cell array containing gene IDs.
%   grRules - cell array containing reaction gene rules.
% Output:
%   rxnGeneMat - logical sparse reaction-gene association matrix.

    % Validate inputs
    numRxns = numel(grRules);
    numGenes = numel(genes);
    
    % Initialize sparse matrix
    rxnGeneMat = sparse(numRxns, numGenes);
    
    % Regular expression pattern to match gene IDs (adaptable for different naming conventions)
    genePattern = '(\<\w+\>)'; % Matches word tokens separated by word boundaries
    
    for i = 1:numRxns
        ruleStr = grRules{i};
        
        % Skip empty rules
        if isempty(ruleStr)
            continue;
        end
        
        % Extract all unique gene IDs
        [geneTokens, ~] = regexp(ruleStr, genePattern, 'tokens', 'match');
        uniqueGenes = unique([geneTokens{:}]);
        
        % Filter out logical operators 'and' and 'or'
        uniqueGenes = uniqueGenes(~ismember(lower(uniqueGenes), {'and','or'}));
        
        % Find the indices of genes in the provided gene list
        [found, geneIdx] = ismember(uniqueGenes, genes);
        
        if any(~found)
            warning('Reaction %d contains unrecognized gene(s): %s', i, strjoin(uniqueGenes(~found), ', '));
            geneIdx = geneIdx(found);
        end
        
        % Update the association matrix
        if ~isempty(geneIdx)
            rxnGeneMat(i, geneIdx) = 1;
        end
    end
    
    % Convert to logical sparse matrix
    % rxnGeneMat = logical(rxnGeneMat);
    
    % Validate matrix dimensions
    assert(size(rxnGeneMat,1) == numRxns, 'Number of reactions does not match.');
    assert(size(rxnGeneMat,2) == numGenes, 'Number of genes does not match.');
end