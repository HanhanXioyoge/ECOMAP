function [ecModel, sluiceConfig] = applySluiceStructure(ecModel, ex_rxn_list, prot_pool)
% applySluiceStructure
%   Applies the sluice structure to exchange reactions BEFORE calibration.
%   This creates a model with dual-branch structure (Basal & Extended).
%
% Input:
%   ecModel      - Base enzyme-constrained model
%   ex_rxn_list - Cell array of exchange reaction names to apply sluice
%                 e.g., {'EX_glc__D_e', 'EX_ac_e', 'EX_succ_e'}
%   prot_pool    - Protein pool name (default: 'prot_pool')
%
% Output:
%   ecModel       - Model with sluice structure applied
%   sluiceConfig  - Structure storing the configuration for later use
%
% Usage:
%   1. Apply structure BEFORE calibration:
%      [model_sluice, config] = applySluiceStructure(model, {'EX_glc__D_e'});
%
%   2. Use in Bayesian calibration with saved kcat sets:
%      model = applySluiceStructure(baseModel, {'EX_glc__D_e'});
%      model.enzymeConstraints.kcat = savedKcats;
%      model = UpdateSmatrix(model);
%
% Note: This function does NOT set Umin values. Umin will be calibrated
%       later by GAUKS or can be set manually using setUmin().

    if nargin < 3 || isempty(prot_pool)
        prot_pool = 'prot_pool';
    end

    if ischar(ex_rxn_list)
        ex_rxn_list = {ex_rxn_list};
    end

    fprintf('[applySluiceStructure] Applying sluice structure to %d reactions...\n', ...
        length(ex_rxn_list));

    % Initialize config storage
    sluiceConfig = struct();
    sluiceConfig.appliedReactions = ex_rxn_list;
    sluiceConfig.umin = zeros(length(ex_rxn_list), 1);  % Placeholder

    % Apply structure to each exchange reaction
    for i = 1:length(ex_rxn_list)
        ex_rxn = ex_rxn_list{i};

        % Check if reaction exists
        rxnIdx = find(strcmp(ecModel.rxns, ex_rxn));
        if isempty(rxnIdx)
            warning('Reaction %s not found in model, skipping.', ex_rxn);
            continue;
        end

        % Find the metabolite (substrate) this exchange consumes
        met_idx = find(ecModel.S(:, rxnIdx) < 0, 1);
        if isempty(met_idx)
            warning('No substrate found for %s, skipping.', ex_rxn);
            continue;
        end

        met_e = ecModel.mets{met_idx};
        met_name = ecModel.metNames{met_idx};

        % Get metabolite properties
        met_form = '';
        if isfield(ecModel, 'metFormulas') && ~isempty(ecModel.metFormulas)
            if met_idx <= length(ecModel.metFormulas)
                met_form = ecModel.metFormulas{met_idx};
            end
        end

        met_charge = 0;
        if isfield(ecModel, 'metCharges')
            if met_idx <= length(ecModel.metCharges)
                met_charge = ecModel.metCharges(met_idx);
            end
        end

        met_comp = 2;
        if isfield(ecModel, 'metComps')
            if met_idx <= length(ecModel.metComps)
                met_comp = ecModel.metComps(met_idx);
            end
        end

        % Create virtual metabolite pool
        met_v = [met_e '_v'];

        % Add virtual metabolite
        ecModel = addMetaboliteEC(ecModel, met_v, met_comp, ...
            'Name', [met_name ' V-Pool'], ...
            'Formula', met_form, 'Charge', met_charge);

        % Store original bounds
        orig_lb = ecModel.lb(rxnIdx);
        orig_ub = ecModel.ub(rxnIdx);

        % Step 1: Redirect original EX reaction to point to virtual pool (Master Switch)
        % Original: met_e -> (cell), Now: met_e -> met_v (virtual pool)
        % Stoichiometry: met_e (outside) is consumed (-1), met_v is produced (+1)
        ecModel = addReactionEC(ecModel, ex_rxn, {met_v}, -1, ...
            orig_lb, orig_ub, 'OnDuplicate', 'overwrite');

        % Step 2: Add Basal Branch - free uptake up to Umin (default: 0, will be calibrated)
        % Reaction: met_v -> met_e (from virtual pool to outside)
        % Stoichiometry: met_v is consumed (-1), met_e is produced (+1)
        ecModel = addReactionEC(ecModel, [ex_rxn '_basal'], ...
            {met_v, met_e}, [1, -1], -1000, 1000);

        % Step 3: Add Extended Branch - penalized uptake beyond Umin
        % Reaction: met_v + prot_pool -> met_e
        % Stoichiometry: met_v consumed (-1), prot_pool consumed (-xi), met_e produced (+1)
        ecModel = addReactionEC(ecModel, [ex_rxn '_extended'], ...
            {met_v, met_e, prot_pool}, [1, -1, 1], -1000, 0);

        fprintf('  Applied: %s -> %s_basal + %s_extended\n', ex_rxn, ex_rxn, ex_rxn);

        % Store in config
        sluiceConfig.umin(i) = 0;  % Placeholder, to be calibrated
    end

    fprintf('[applySluiceStructure] Done. Model now has sluice structure.\n');
end
