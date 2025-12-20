function [model, tunedKcats] = sensitivityTuning(model, desiredGrowthRate, foldChange, protToIgnore, verbose, parameters)
% sensitivityTuning
%    Function that relaxes the most limiting kcats until a certain growth rate
%    is reached. The function will update kcats in model.enzymeConstraints.kcat.
%
% Input:
%   model              an ecModel in ECOMAP format
%   desiredGrowthRate  kcats will be relaxed until this growth rate is reached
%   parameters         a structure containing model parameters, including
%                      bioRxn, prot_pool, and other model-specific info
%   foldChange         kcat values will be increased by this fold-change.
%                      (Opt, default 10)
%   protToIgnore       vector of protein ids to be ignored in tuned kcats.
%                      e.g. {'P38122', 'Q99271'} (Optional, default = [])
%   verbose            logical whether progress should be reported (Optional,
%                      default true)
%
% Output:
%   model              ecModel with updated model.enzymeConstraints.kcat
%   tunedKcats         structure with information on tuned kcat values
%                      rxns     identifiers of reactions with tuned kcat
%                               values
%                      rxnNames names of the reactions in tunedKcats.rxns
%                      enzymes  enzymes that catalyze the reactions in
%                               tunedKcats.rxns, whose kcat value has been
%                               tuned.
%                      oldKcat  kcat values in the input model
%                      newKcat  kcat values in the output model, after tuning
%
% From GECKO:  https://github.com/SysBioChalmers/GECKO/blob/main/src/geckomat/kcat_sensitivity_analysis/sensitivityTuning.m


% Default input arguments
if nargin < 6 || isempty(parameters)
    parameters = ParameterManager.getParams();
    if isempty(parameters), error('ParameterManager is not set.'); end
end

if nargin < 5 || isempty(verbose)
    verbose = true;
end
if nargin < 4 || isempty(protToIgnore)
    protToIgnore = {};
end
if nargin < 3 || isempty(foldChange)
    foldChange = 10;
end

% Initialize model and parameters based on type (basic, isozyme, integrated)
kcatList = [];
m = model;

% Ensure growth maximization
m.c = double(strcmp(m.rxns, parameters.bioRxn));

% Solve linear program for growth rate
[res, hs] = solveLP(m);
if isempty(res.x)
    error('FBA of input model gives no valid result. Check protein pool constraints and exchange constraints.')
end

lastGrowth = 0;

% Handle different model types (basic, isozyme, integrated)
if strcmp(model.enzymeConstraints.ecModeltype, 'integrated')
    % Full model type - process by reaction draw and protein pool
    drawRxns = startsWith(m.rxns, 'usage_prot_');
    idxToIgnore = cellfun(@(x) find(strcmpi(model.rxns, ['usage_prot_' x])), protToIgnore);
    iteration = 1;
    
    while true
        [res, hs] = solveLP(m, 0, [], hs); % Skip parsimonius, just time-consuming
        if (lastGrowth == res.f)
            printOrange('WARNING: No growth increase from increased kcats - check uptake reaction constraints.\n');
            break;
        end
        lastGrowth = res.f;
        if verbose; disp(['Iteration ' num2str(iteration) ': Growth: ' num2str(lastGrowth)]); end
        if (lastGrowth >= desiredGrowthRate)
            break;
        end
        
        iteration = iteration + 1;
        
        % Find the highest draw_prot rxn flux
        drawFluxes = zeros(length(drawRxns), 1);
        drawFluxes(drawRxns) = res.x(drawRxns);
        
        % Remove proteins to ignore
        drawFluxes(idxToIgnore) = 0;
        [x, sel] = min(drawFluxes); % Select minimum flux
        metSel = m.S(:, sel) < 0;   % Negative coefficients (consumption)
        
        % Find the reaction with the largest protein consumption
        protFluxes = m.S(metSel, :).' .* res.x;
        [~, rxnSel] = min(protFluxes);
        
        % Update kcat values for selected reactions
        kcatList = [kcatList, rxnSel];
        rxn = m.rxns(rxnSel);
        targetSubRxn = strcmp(m.enzymeConstraints.rxns, rxn);
        
        % Update notes for the reaction being tuned
        if ~strcmp(m.enzymeConstraints.source(targetSubRxn), 'sensitivityTuning')
            oldNote = m.enzymeConstraints.notes{targetSubRxn};
            newNote = ['preTuneKcat=' num2str(m.enzymeConstraints.kcat(targetSubRxn)) ' | source:' m.enzymeConstraints.source{targetSubRxn}];
            if ~isempty(oldNote) || ~ strcmp(oldNote, 'sensitivityTuning')
                newNote = [oldNote '; ' newNote];
            end
            m.enzymeConstraints.notes{targetSubRxn} = newNote;
        end
        m.enzymeConstraints.kcat(targetSubRxn) = m.enzymeConstraints.kcat(targetSubRxn) * foldChange;
        m.enzymeConstraints.source(targetSubRxn) = {'sensitivityTuning'};
        m = UpdateSmatrix(m, targetSubRxn);
    end
elseif strcmp(model.enzymeConstraints.ecModeltype, 'isozyme')
    origRxns = m.enzymeConstraints.rxns;                  % For isozyme models
    idxToIgnore = cellfun(@(x) find(m.enzymeConstraints.rxnEnzMat(:, strcmpi(m.enzymeConstraints.enzymes, x))), protToIgnore, 'UniformOutput', false);
    idxToIgnore = unique(cat(1, idxToIgnore{:}));
    idxToIgnore = cellfun(@(x) find(strcmpi(m.rxns, x)), origRxns(idxToIgnore));
    iteration = 1;
    while true
        res = solveLP(m, 0); % Skip parsimonius, only time-consuming
        if (lastGrowth == res.f)
            printOrange('No growth increase from increased kcats - check uptake reaction constraints.\n');
            break;
        end
        lastGrowth = res.f;
        if verbose; disp(['Iteration ' num2str(iteration) ': Growth: ' num2str(lastGrowth)]); end
        if (lastGrowth >= desiredGrowthRate)
            break;
        end
        
        iteration = iteration + 1;
        
        % Find the highest protein usage flux
        protPoolStoich = m.S(strcmp(m.mets, 'prot_pool'), :).' ;
        protPoolStoich(idxToIgnore) = 0;
        [~, sel] = min(res.x .* protPoolStoich);  % Max consumption
        
        kcatList = [kcatList, sel];
        rxn = m.rxns(sel.');
        targetSubRxn = strcmp(origRxns, rxn);
        % Update notes for the reaction being tuned

        if ~strcmp(m.enzymeConstraints.source(targetSubRxn), 'sensitivityTuning')
            oldNote = m.enzymeConstraints.notes{targetSubRxn};
            newNote = ['preTuneKcat=' num2str(m.enzymeConstraints.kcat(targetSubRxn)) ' | source:' m.enzymeConstraints.source{targetSubRxn}];
            if ~isempty(oldNote) || ~ strcmp(oldNote, 'sensitivityTuning')
                newNote = [oldNote '; ' newNote];
            end
            m.enzymeConstraints.notes{targetSubRxn} = newNote;
        end
        % Update kcat values for selected reactions
        m.enzymeConstraints.kcat(targetSubRxn) = m.enzymeConstraints.kcat(targetSubRxn) * foldChange;
        m = UpdateSmatrix(m, rxn);
    end
elseif strcmp(model.enzymeConstraints.ecModeltype, 'basic')
    origRxns = extractAfter(m.enzymeConstraints.rxns, 4); % For basic models
    idxToIgnore = cellfun(@(x) find(m.enzymeConstraints.rxnEnzMat(:, strcmpi(m.enzymeConstraints.enzymes, x))), protToIgnore, 'UniformOutput', false);
    idxToIgnore = unique(cat(1, idxToIgnore{:}));
    idxToIgnore = cellfun(@(x) find(strcmpi(m.rxns, x)), origRxns(idxToIgnore));
    
    iteration = 1;
    while true
        res = solveLP(m, 0); % Skip parsimonius, only time-consuming
        if (lastGrowth == res.f)
            printOrange('No growth increase from increased kcats - check uptake reaction constraints.\n');
            break;
        end
        lastGrowth = res.f;
        if verbose; disp(['Iteration ' num2str(iteration) ': Growth: ' num2str(lastGrowth)]); end
        if (lastGrowth >= desiredGrowthRate)
            break;
        end

        iteration = iteration + 1;
        
        % Find the highest protein usage flux
        protPoolStoich = m.S(strcmp(m.mets, 'prot_pool'), :).' ;
        protPoolStoich(idxToIgnore) = 0;
        [~, sel] = min(res.x .* protPoolStoich);  % Max consumption
        
        kcatList = [kcatList, sel];
        rxn = m.rxns(sel.');
        targetSubRxns = strcmp(origRxns, rxn);
        % Update notes for the reaction being tuned

        if ~strcmp(m.enzymeConstraints.source(targetSubRxns), 'sensitivityTuning')
            m.enzymeConstraints.notes(targetSubRxns) = {'sensitivityTuning'};
        end
        % Update kcat values for selected reactions
        m.enzymeConstraints.kcat(targetSubRxns) = m.enzymeConstraints.kcat(targetSubRxns) * foldChange;
        m = UpdateSmatrix(m, rxn);
    end
end

% Collect results in tunedKcats
kcatList = unique(kcatList);
tunedKcats.rxns = m.rxns(kcatList);
tunedKcats.rxnNames = m.rxnNames(kcatList);

% Find enzymes corresponding to the reactions with tuned kcats
if strcmp(model.enzymeConstraints.ecModeltype, 'integrated')
    [~, rxnIdx] = ismember(tunedKcats.rxns, m.enzymeConstraints.rxns);
else
    [~, rxnIdx] = ismember(tunedKcats.rxns, origRxns);
end

tunedKcats.enzymes = cell(numel(kcatList), 1);
for i = 1:numel(rxnIdx)
    [~, metIdx] = find(m.enzymeConstraints.rxnEnzMat(rxnIdx(i), :));
    tunedKcats.enzymes{i} = strjoin(m.enzymeConstraints.enzymes(metIdx), ';');
end
tunedKcats.oldKcat = model.enzymeConstraints.kcat(rxnIdx);
tunedKcats.newKcat = m.enzymeConstraints.kcat(rxnIdx);
tunedKcats.source = model.enzymeConstraints.source(rxnIdx);

model = m;

end