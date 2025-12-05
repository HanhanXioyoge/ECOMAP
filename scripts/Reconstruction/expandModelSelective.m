function [newModel, rxnToCheck] = expandModelSelective(model, expandMaskOrIdx)
% expandModelSelective
%   Expands ONLY the reactions selected by expandMaskOrIdx when their GPR contains ' or '.
%   Other reactions are kept untouched even if they have ' or ' in grRules.
%
% Inputs:
%   model            : COBRA-like model structure
%   expandMaskOrIdx  : logical [nRxn x 1] or vector of indices to be eligible for expansion
%
% Outputs:
%   newModel         : model with selected reactions split by isozyme (' or ') into copies
%                      with IDs suffixed as _EXP_1, _EXP_2, ...
%   rxnToCheck       : cellstr of original reaction IDs that contained nested 'and'/'or'
%                      relationships (recommended for manual inspection)
%
% Notes:
%   - This is a drop-in selective version of your original expandModel; the only change is
%     that we zero out numOrs for reactions not in the selection, so they never expand.

    if nargin < 2 || isempty(expandMaskOrIdx)
        % default: expand all candidates (behaves like the original expandModel)
        expandMask = true(numel(model.rxns),1);
    else
        if islogical(expandMaskOrIdx)
            expandMask = expandMaskOrIdx(:);
        else
            expandMask = false(numel(model.rxns),1);
            expandMask(expandMaskOrIdx(:)) = true;
        end
    end

    % Pre-count ' or ' per reaction, then ZERO OUT for non-selected reactions
    numOrs = count(model.grRules, ' or ');
    numOrs(~expandMask) = 0;

    toAdd       = sum(numOrs);
    prevNumRxns = numel(model.rxns);
    rxnToCheck  = {};

    if toAdd > 0
        % Indices to copy (each reaction with k "or" adds k copies)
        cpyIndices = repelem(1:prevNumRxns, numOrs);

        % Copy simple-aligned fields
        model.S         = [model.S, model.S(:, cpyIndices)];
        model.rxnNames  = [model.rxnNames;  model.rxnNames(cpyIndices)];
        model.lb        = [model.lb;        model.lb(cpyIndices)];
        model.ub        = [model.ub;        model.ub(cpyIndices)];
        model.rev       = [model.rev;       model.rev(cpyIndices)];
        model.c         = [model.c;         model.c(cpyIndices)];
        if isfield(model,'subSystems')
            model.subSystems = [model.subSystems; model.subSystems(cpyIndices)];
        end
        if isfield(model,'eccodes')
            model.eccodes = [model.eccodes; model.eccodes(cpyIndices)];
        end
        if isfield(model,'equations')
            model.equations = [model.equations; model.equations(cpyIndices)];
        end
        if isfield(model,'rxnMiriams')
            model.rxnMiriams = [model.rxnMiriams; model.rxnMiriams(cpyIndices)];
        end
        if isfield(model,'rxnComps')
            model.rxnComps = [model.rxnComps; model.rxnComps(cpyIndices)];
        end
        if isfield(model,'rxnFrom')
            model.rxnFrom = [model.rxnFrom; model.rxnFrom(cpyIndices)];
        end
        if isfield(model,'rxnNotes')
            model.rxnNotes = [model.rxnNotes; model.rxnNotes(cpyIndices)];
        end
        if isfield(model,'rxnReferences')
            model.rxnReferences = [model.rxnReferences; model.rxnReferences(cpyIndices)];
        end
        if isfield(model,'rxnConfidenceScores')
            model.rxnConfidenceScores = [model.rxnConfidenceScores; model.rxnConfidenceScores(cpyIndices)];
        end
        if isfield(model,'rxnDeltaG')
            model.rxnDeltaG = [model.rxnDeltaG; model.rxnDeltaG(cpyIndices)];
        end

        % Prepare variable-sized fields to be filled later
        model.rxns       = [model.rxns;     cell(toAdd,1)];
        model.grRules    = [model.grRules;  cell(toAdd,1)];
        model.rxnGeneMat = [model.rxnGeneMat; sparse(toAdd, size(model.rxnGeneMat,2))];

        % Fill expanded data
        nextIndex = prevNumRxns + 1;
        for i = 1:prevNumRxns
            if numOrs(i) > 0
                % Warn if nested 'and'/'or'
                if ~isempty(strfind(model.grRules{i}, ' and '))
                    rxnToCheck{end+1,1} = model.rxns{i}; %#ok<AGROW>
                end

                geneString = model.grRules{i};
                geneString = strrep(geneString, '(', '');
                geneString = strrep(geneString, ')', '');
                geneString = strrep(geneString, ' or ', ';');

                geneNames = regexp(geneString, ';', 'split');

                % Update reaction i to only use the first gene/unit
                model.grRules{i}  = ['(' geneNames{1} ')'];
                model.rxnGeneMat(i,:) = 0;

                if ~isempty(strfind(geneNames(1), ' and '))
                    andGenes = regexp(geneNames{1}, ' and ', 'split');
                    model.rxnGeneMat(i, ismember(model.genes, andGenes)) = 1;
                else
                    [~, index] = ismember(geneNames(1), model.genes);
                    model.rxnGeneMat(i, index) = 1;
                end

                % Insert additional copies
                for j = 2:numel(geneNames)
                    ind = nextIndex + j - 2;
                    model.rxns{ind} = [model.rxns{i} '_EXP_' num2str(j)];
                    model.grRules{ind} = ['(' geneNames{j} ')'];

                    if ~isempty(strfind(geneNames(j), ' and '))
                        andGenes = regexp(geneNames{j}, ' and ', 'split');
                        model.rxnGeneMat(ind, ismember(model.genes, andGenes)) = 1;
                    else
                        model.rxnGeneMat(ind, ismember(model.genes, geneNames(j))) = 1;
                    end
                end
                model.rxns{i} = [model.rxns{i}, '_EXP_1'];
                nextIndex = nextIndex + numOrs(i);
            end
        end
        newModel = model;
    else
        % Nothing to expand
        newModel = model;
    end

    % Standardize grRules and rebuild rxnGeneMat
    [grRules, rxnGeneMat] = standardizeGrRules(newModel, true);
    newModel.grRules      = grRules;
    newModel.rxnGeneMat   = rxnGeneMat;
end
