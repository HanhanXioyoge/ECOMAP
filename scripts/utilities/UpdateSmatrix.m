function model = UpdateSmatrix(model, updateRxns)
% UpdateSmatrix
%   This function updates the stoichiometric matrix (S) of an ECOMAP ecModel 
%   by applying enzyme constraints based on the provided kcat values. 
%   Depending on the model type (integrated, isozyme, or basic), 
%   it will update the reaction constraints accordingly, considering enzyme 
%   efficiencies, enzyme subunit copies, and molecular weights.
% 
% Inputs:
%   model        - ECOMAP ecModel with enzymeConstraints structure
%   updateRxns   - A logical vector specifying which reactions to update 
%                  (optional, if not provided, all reactions will be updated).
%                  Alternatively, can be a numeric vector, cell array of strings, or 
%                  a string with the reaction identifiers to update.
%
% Outputs:
%   model        - The updated ECOMAP ecModel with the modified stoichiometric matrix S.

% Check if reaction list should be based on enzymeConstraints or model.rxns
if strcmp(model.enzymeConstraints.ecModeltype, 'integrated')
    rxns = model.enzymeConstraints.rxns;  % For 'integrated' model, use enzymeConstraints.rxns
else
    rxns = model.rxns;  % For 'isozyme' and 'basic' models, use model.rxns
end

% Validate updateRxns input and ensure its length matches the number of reactions
if nargin < 2
    updateRxns = true(numel(rxns), 1);  % Default to updating all reactions
elseif isnumeric(updateRxns)
    updateRxnsLog = false(numel(rxns), 1);
    updateRxnsLog(updateRxns) = true;
    updateRxns = updateRxnsLog;  % Convert numeric indices to logical
elseif iscellstr(updateRxns) || ischar(updateRxns) || isstring(updateRxns)
    updateRxnsIds = convertCharArray(updateRxns);
    updateRxns = ismember(rxns, updateRxnsIds);  % Convert reaction IDs to logical vector
end

% Ensure updateRxns is not empty or invalid
if isempty(find(updateRxns, 1)) || isempty(updateRxns)
     error('No reaction to update or updateRxns is logical but without any true value');
end

% Check if the model has enzymeConstraints structure
if ~isfield(model, 'enzymeConstraints')
    error('No model.enzymeConstraints structure could be found: the provided model is not a valid ECOMAP ecModel.');
end

% If no kcat values are provided, exit the function
if all(model.enzymeConstraints.kcat == 0) && isempty(model.enzymeConstraints.kcat)
    warning('WARNING: No kcat values are provided in model.enzymeConstraints.kcat, model remains unchanged.\n');
    return;
end

% Clear existing enzyme usage incorporation in the stoichiometric matrix for 'integrated' model
if strcmp(model.enzymeConstraints.ecModeltype, 'integrated')
    protMetIdx = startsWith(model.mets, 'prot_') & ~strcmp(model.mets, 'prot_pool');
    metabolRxn = unique(model.enzymeConstraints.rxns(updateRxns));
    metabolRxn = ismember(model.rxns, metabolRxn);
    model.S(protMetIdx, metabolRxn) = 0;  % Reset enzyme usage in the matrix
end

% Update kcat values in the stoichiometric matrix depending on the model type
if strcmp(model.enzymeConstraints.ecModeltype, 'integrated')
    % For 'integrated' model: Update reaction constraints considering enzyme subunits
    newKcats = zeros(numel(updateRxns) * 10, 5);  % Pre-allocate for kcat data
    updateRxns = find(updateRxns);  % Find reactions to update
    kcatFirst = 0;
    
    for i = 1:numel(updateRxns)
        j = updateRxns(i);
        if model.enzymeConstraints.kcat(j) ~= 0
            enzymes = find(model.enzymeConstraints.rxnEnzMat(j, :));  % Get enzyme indices for reaction
            kcatLast = kcatFirst + numel(enzymes);
            kcatFirst = kcatFirst + 1;
            newKcats(kcatFirst:kcatLast, 1) = j;
            newKcats(kcatFirst:kcatLast, 2) = enzymes;
            newKcats(kcatFirst:kcatLast, 3) = model.enzymeConstraints.rxnEnzMat(j, enzymes);
            newKcats(kcatFirst:kcatLast, 4) = model.enzymeConstraints.kcat(j);
            newKcats(kcatFirst:kcatLast, 5) = model.enzymeConstraints.mw(enzymes);
            kcatFirst = kcatLast;
        end
    end
    
    if exist('kcatLast', 'var')
        newKcats(kcatLast + 1:end, :) = [];  % Remove empty rows

        sel = newKcats(:, 4) > 0;  % Only apply non-zero kcat values
        newKcats(sel, 4) = newKcats(sel, 4) * 3600;  % Convert kcat from per second to per hour
        newKcats(sel, 4) = newKcats(sel, 5) ./ newKcats(sel, 4);  % Calculate MW/kcat ratio
        newKcats(sel, 4) = newKcats(sel, 3) .* newKcats(sel, 4);  % Multiply by subunit copies
        newKcats(~sel, 4) = 0;  % Set zero for invalid kcat values

        % Update the stoichiometric matrix with enzyme usage
        [~, newKcats(:, 1)] = ismember(model.enzymeConstraints.rxns(newKcats(:, 1)), model.rxns);
        [~, newKcats(:, 2)] = ismember(strcat('prot_', model.enzymeConstraints.enzymes(newKcats(:, 2))), model.mets);
        linearIndices = sub2ind(size(model.S), newKcats(:, 2), newKcats(:, 1));
        model.S(linearIndices) = -newKcats(:, 4);  % Update the stoichiometric matrix with the calculated kcat values
    end

elseif strcmp(model.enzymeConstraints.ecModeltype, 'isozyme')
    % For 'isozyme' model: Update enzyme usage based on MW/kcat values for the enzyme(s) involved in the reaction
    prot_pool_idx = find(ismember(model.mets, 'prot_pool'));
    modRxns = model.enzymeConstraints.rxns;
    [hasEc, ~] = ismember(model.rxns, modRxns);
    hasEc = find(hasEc & updateRxns);  % Get reactions to update
    [~, rxnIdx] = ismember(modRxns, model.rxns);
    
    for i = 1:numel(hasEc)
        % For each reaction, get the enzyme(s) involved (one enzyme per reaction in isozyme model)
        ecIdx = find(rxnIdx == hasEc(i));
        % Calculate MW/kcat for the enzyme(s) involved in the reaction
        MWkcat = (model.enzymeConstraints.rxnEnzMat(ecIdx, :) * model.enzymeConstraints.mw) ./ model.enzymeConstraints.kcat(ecIdx);
        MWkcat(isinf(MWkcat) | isnan(MWkcat)) = 0;  % Handle invalid MW/kcat values
        % Update the stoichiometric matrix for the reaction based on MW/kcat
        model.S(prot_pool_idx, hasEc(i)) = -sum(MWkcat / 3600);  % Convert to per hour and update the stoichiometry
    end

elseif strcmp(model.enzymeConstraints.ecModeltype, 'basic')
    % For 'basic' model: Update enzyme usage based on MW/kcat, using the most efficient enzyme
    prot_pool_idx = find(ismember(model.mets, 'prot_pool'));
    modRxns = extractAfter(model.enzymeConstraints.rxns, 4);
    [hasEc, ~] = ismember(model.rxns, modRxns);
    hasEc = find(hasEc & updateRxns);  % Get reactions to update
    [~, rxnIdx] = ismember(modRxns, model.rxns);
    
    for i = 1:numel(hasEc)
        ecIdx = find(rxnIdx == hasEc(i));
        MWkcat = (model.enzymeConstraints.rxnEnzMat(ecIdx, :) * model.enzymeConstraints.mw) ./ model.enzymeConstraints.kcat(ecIdx);
        MWkcat(isinf(MWkcat) | isnan(MWkcat)) = 0;  % Handle invalid MW/kcat values
        model.S(prot_pool_idx, hasEc(i)) = -min(MWkcat / 3600);  % Select the most efficient enzyme and update the stoichiometric matrix
    end
end

end
