function [model, tunedKcats] = sensitivityTuning(model, desiredGrowthRate, foldChange, protToIgnore, verbose, parameters)
% sensitivityTuning
%    Relaxes the most limiting kcat values until a target growth rate is reached.
%    This version uses a parsimonious step to ensure only true bottlenecks are tuned.

% 1. Default argument handling
if nargin < 6 || isempty(parameters)
    parameters = ParameterManager.getParams();
end
if nargin < 5 || isempty(verbose), verbose = true; end
if nargin < 4 || isempty(protToIgnore), protToIgnore = {}; end
if nargin < 3 || isempty(foldChange), foldChange = 10; end

kcatList = [];
m = model;
m.c = double(strcmp(m.rxns, parameters.bioRxn)); % Set growth as objective

% Initial FBA check
[res, hs] = solveLP(m);
if isempty(res.x)
    error('Initial FBA failed. Check protein pool or exchange constraints.');
end

lastGrowth = 0;
iteration = 1;

% 2. Process based on ecModel type
if strcmp(model.enzymeConstraints.ecModeltype, 'integrated')
    % Logic for Integrated Models (usage_prot_ reactions)
    drawRxns = startsWith(m.rxns, 'usage_prot_');
    idxToIgnore = cellfun(@(x) find(strcmpi(model.rxns, ['usage_prot_' x])), protToIgnore);
    
    while true
        % --- Step A: Maximize Growth ---
        res = solveLP(m, 0, [], hs); 
        if isempty(res.x) || (lastGrowth == res.f && iteration > 1)
            if verbose; disp('No further growth increase possible.'); end
            break;
        end
        lastGrowth = res.f;
        if verbose; fprintf('Iteration %d: Growth = %.4f\n', iteration, lastGrowth); end
        if (lastGrowth >= desiredGrowthRate), break; end
        
        % --- Step B: Parsimonious Step (Minimize Protein Pool) ---
        % Fix growth rate to identify the most efficient protein allocation
        m_pfba = m;
        idxBio = find(strcmp(m_pfba.rxns, parameters.bioRxn));
        m_pfba.lb(idxBio) = lastGrowth * 0.999;
        m_pfba.ub(idxBio) = lastGrowth * 1.001;
        m_pfba.c(:) = 0;
        protExIdx = strcmp(m_pfba.rxns, 'prot_pool_exchange');
        m_pfba.c(protExIdx) = 1; % Maximize exchange (effectively minimizing consumption)
        
        res_pfba = solveLP(m_pfba);
        if isempty(res_pfba.x), res_final = res; else, res_final = res_pfba; end
        
        % --- Step C: Identify and Tune Bottleneck ---
        drawFluxes = zeros(length(drawRxns), 1);
        drawFluxes(drawRxns) = res_final.x(drawRxns);
        drawFluxes(idxToIgnore) = 0;
        
        [~, sel] = min(drawFluxes); % Highest consumption (most negative flux)
        metSel = m.S(:, sel) < 0;   % Find the enzyme metabolite
        
        % Identify the specific enzymatic reaction causing the draw
        protFluxes = m.S(metSel, :).' .* res_final.x;
        [~, rxnSel] = min(protFluxes);
        
        % Update kcat
        kcatList = [kcatList, rxnSel];
        targetIdx = strcmp(m.enzymeConstraints.rxns, m.rxns(rxnSel));
        m = updateKcat(m, targetIdx, foldChange);
        iteration = iteration + 1;
    end

elseif any(strcmp(model.enzymeConstraints.ecModeltype, {'isozyme', 'basic'}))
    % Logic for Isozyme or Basic Models
    isIso = strcmp(model.enzymeConstraints.ecModeltype, 'isozyme');
    if isIso
        origRxns = m.enzymeConstraints.rxns;
    else
        origRxns = extractAfter(m.enzymeConstraints.rxns, 4);
    end
    
    % Map proteins to ignore to reaction indices
    idxToIgnore = [];
    for i = 1:numel(protToIgnore)
        enzIdx = find(strcmpi(m.enzymeConstraints.enzymes, protToIgnore{i}));
        if ~isempty(enzIdx)
            [rxnRows, ~] = find(m.enzymeConstraints.rxnEnzMat(:, enzIdx));
            for j = 1:numel(rxnRows)
                idxToIgnore(end+1) = find(strcmpi(m.rxns, origRxns{rxnRows(j)}));
            end
        end
    end

    while true
        % --- Step A: Maximize Growth ---
        res = solveLP(m, 0);
        if isempty(res.x) || (lastGrowth == res.f && iteration > 1), break; end
        lastGrowth = res.f;
        if verbose; fprintf('Iteration %d: Growth = %.4f\n', iteration, lastGrowth); end
        if (lastGrowth >= desiredGrowthRate), break; end

        % --- Step B: Parsimonious Step ---
        m_pfba = m;
        idxBio = find(strcmp(m_pfba.rxns, parameters.bioRxn));
        m_pfba.lb(idxBio) = lastGrowth * 0.999;
        m_pfba.ub(idxBio) = lastGrowth * 1.001;
        m_pfba.c(:) = 0;
        m_pfba.c(strcmp(m_pfba.rxns, 'prot_pool_exchange')) = 1;
        res_pfba = solveLP(m_pfba);
        if isempty(res_pfba.x), res_final = res; else, res_final = res_pfba; end

        % --- Step C: Identify Bottleneck ---
        protPoolStoich = m.S(strcmp(m.mets, 'prot_pool'), :).';
        protPoolStoich(idxToIgnore) = 0;
        [~, sel] = min(res_final.x .* protPoolStoich); 
        
        kcatList = [kcatList, sel];
        rxnID = m.rxns{sel};
        targetIdx = strcmp(m.enzymeConstraints.rxns, rxnID);
        if ~any(targetIdx) && ~isIso % Basic model prefix handling
            targetIdx = strcmp(m.enzymeConstraints.rxns, ['enz_' rxnID]);
        end
        
        m = updateKcat(m, targetIdx, foldChange);
        iteration = iteration + 1;
    end
end

% 3. Format Results
kcatList = unique(kcatList);
tunedKcats.rxns = m.rxns(kcatList);
tunedKcats.rxnNames = m.rxnNames(kcatList);

if strcmp(model.enzymeConstraints.ecModeltype, 'integrated')
    [~, rxnIdxInEC] = ismember(tunedKcats.rxns, m.enzymeConstraints.rxns);
else
    [~, rxnIdxInEC] = ismember(tunedKcats.rxns, origRxns);
end

% Map enzymes for the tuned reactions
tunedKcats.enzymes = cell(numel(kcatList), 1);
for i = 1:numel(rxnIdxInEC)
    [~, enzIdxs] = find(m.enzymeConstraints.rxnEnzMat(rxnIdxInEC(i), :));
    tunedKcats.enzymes{i} = strjoin(m.enzymeConstraints.enzymes(enzIdxs), '; ');
end

tunedKcats.oldKcat = model.enzymeConstraints.kcat(rxnIdxInEC);
tunedKcats.newKcat = m.enzymeConstraints.kcat(rxnIdxInEC);
tunedKcats.source = m.enzymeConstraints.source(rxnIdxInEC);
model = m;

end

%% --- Internal Helper: Update kcat and Metadata ---
function m = updateKcat(m, targetIdx, foldChange)
    if ~strcmp(m.enzymeConstraints.source(targetIdx), 'sensitivityTuning')
        oldNote = m.enzymeConstraints.notes{targetIdx};
        newNote = sprintf('preTuneKcat=%g | source:%s', ...
            m.enzymeConstraints.kcat(targetIdx), m.enzymeConstraints.source{targetIdx});
        if ~isempty(oldNote)
            m.enzymeConstraints.notes{targetIdx} = [oldNote '; ' newNote];
        else
            m.enzymeConstraints.notes{targetIdx} = newNote;
        end
    end
    m.enzymeConstraints.kcat(targetIdx) = m.enzymeConstraints.kcat(targetIdx) * foldChange;
    
    % Efficiently update the stoichiometric matrix for the modified kcat
    rxnID = m.enzymeConstraints.rxns{targetIdx};
    m = UpdateSmatrix(m, rxnID);
end