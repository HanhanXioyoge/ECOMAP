function reducedModel = removemodelGenes(model, genesToRemove, standardizeRules)
% removeGenes
%   Deletes a set of genes from a model and modifies related gene rules (gpr),
%   gene names, gene Miriams, and gene short names. It keeps the reactions
%   and flux bounds unchanged.
%
%   model                   a model structure
%   genesToRemove           a cell array of gene IDs to remove
%   standardizeRules        format gene rules to be compliant with standard format
%
%   reducedModel            an updated model structure

if nargin < 3
    standardizeRules = true;
end

% Format grRules and rxnGeneMatrix:
if standardizeRules
    [grRules, rxnGeneMat, ~] = standardizeGrRules(model, true);
    model.grRules = grRules;
    model.rxnGeneMat = rxnGeneMat;
else
    rxnGeneMat = model.rxnGeneMat;
end

reducedModel = model;

% Only remove genes that are actually in the model
if ~(islogical(genesToRemove) || isnumeric(genesToRemove))
    genesToRemove = convertCharArray(genesToRemove);
    genesToRemove = genesToRemove(ismember(genesToRemove, model.genes));
end

if ~isempty(genesToRemove)
    indexesToRemove = getIndexes(model, genesToRemove, 'genes');
    
    if ~isempty(indexesToRemove)
        % Make corresponding columns in rxnGeneMat 0:
        reducedModel.rxnGeneMat(:, indexesToRemove) = 0;        
        genes = model.genes(indexesToRemove);

        % Update grRules for each reaction that contains the gene to be removed
        for i = 1:length(genes)
            % Find all reactions for this gene and loop through them:
            geneRxns = find(rxnGeneMat(:, indexesToRemove(i)));
            if ~isempty(geneRxns)
                for j = 1:numel(geneRxns)
                    index = geneRxns(j);
                    grRule = reducedModel.grRules{index};
                    if ~isempty(grRule)
                        % Adapt gene rule & gene matrix:
                        grRule = removeGeneFromRule(grRule, genes{i});
                        reducedModel.grRules{index} = grRule;
                    end
                end
            end
        end
        
        % Update genes, geneMiriams, and geneShortNames
        reducedModel.genes(indexesToRemove) = [];
        reducedModel.geneMiriams(indexesToRemove) = [];
        reducedModel.geneShortNames(indexesToRemove) = [];
    end
end

% Format grRules and rxnGeneMatrix after all modifications
if standardizeRules
    [grRules, rxnGeneMat] = standardizeGrRules(reducedModel, true);
    reducedModel.grRules = grRules;
    reducedModel.rxnGeneMat = rxnGeneMat;
end

end

function geneRule = removeGeneFromRule(geneRule, geneToRemove)
% This function receives a standard gene rule and it returns it without the
% chosen gene.
geneSets = strsplit(geneRule, ' or ');
hasGene = ~cellfun(@isempty, strfind(geneSets, geneToRemove));
geneSets = geneSets(~hasGene);
geneRule = strjoin(geneSets, ' or ');
end
