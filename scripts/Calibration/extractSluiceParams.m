function sluiceParams = extractSluiceParams(ecModel, ex_rxn_list, prot_pool)
% EXTRACTSLUICEPARAMS
%   Extract Umin and Xi values from a model with sluice structure.
%
% Input:
%   ecModel     - Model with sluice structure (from applySluiceStructure)
%   ex_rxn_list - Cell array of exchange reaction names
%   prot_pool   - Protein pool name (default: 'prot_pool')
%
% Output:
%   sluiceParams - struct with fields {reactions, umin, xi}
%
% Usage:
%   % Extract parameters from calibrated models
%   ex_rxn_list = {'EX_glc__D_e', 'EX_ac_e', ...};
%   sluiceParams = extractSluiceParams(calibratedModel, ex_rxn_list);
%
%   % Save to KcatRepo
%   kcatRepo.addGroup('Bayesian_FullData', kcats, sluiceParams, 'Bayesian', 'dataTag');

    if nargin < 3 || isempty(prot_pool)
        prot_pool = 'prot_pool';
    end

    if ischar(ex_rxn_list)
        ex_rxn_list = {ex_rxn_list};
    end

    nRxns = length(ex_rxn_list);
    umin = zeros(nRxns, 1);
    xi = zeros(nRxns, 1);

    fprintf('[extractSluiceParams] Extracting %d reactions...\n', nRxns);

    for i = 1:nRxns
        ex_rxn = ex_rxn_list{i};

        % Get Umin from basal branch bounds
        basal_rxn = [ex_rxn '_basal'];
        basalIdx = find(strcmp(ecModel.rxns, basal_rxn));

        if isempty(basalIdx)
            warning('Basal reaction %s not found', basal_rxn);
            umin(i) = NaN;
            xi(i) = NaN;
            continue;
        end

        % Umin is the negative of lb (uptake rate is negative)
        umin(i) = abs(ecModel.lb(basalIdx));

        % Get Xi from extended branch S matrix
        ext_rxn = [ex_rxn '_extended'];
        extIdx = find(strcmp(ecModel.rxns, ext_rxn));

        if isempty(extIdx)
            warning('Extended reaction %s not found', ext_rxn);
            xi(i) = NaN;
            continue;
        end

        protRow = find(strcmp(ecModel.mets, prot_pool));
        if ~isempty(protRow)
            xi(i) = ecModel.S(protRow, extIdx);
        else
            xi(i) = NaN;
        end

        fprintf('  %s: Umin=%.4f, Xi=%.4f\n', ex_rxn, umin(i), xi(i));
    end

    % Create struct with proper dimensions
    sluiceParams = struct();
    sluiceParams.reactions = ex_rxn_list(:)';  % Ensure row vector
    sluiceParams.umin = umin(:)';               % Ensure row vector
    sluiceParams.xi = xi(:)';

    fprintf('[extractSluiceParams] Done.\n');
end
