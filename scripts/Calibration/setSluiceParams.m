function ecModel = setSluiceParams(ecModel, ex_rxn, umin, xi, prot_pool)
% SETSLUICEPARAMS
%   Sets the Umin and Xi (proteomic cost) values for a sluice-structured exchange reaction.
%
% Input:
%   ecModel     - Model with sluice structure (from applySluiceStructure)
%   ex_rxn     - Exchange reaction name (e.g., 'EX_glc__D_e')
%   umin       - Minimum uptake rate (umol/gDW/h)
%   xi         - Proteomic cost coefficient (for extended branch)
%   prot_pool   - Protein pool name (default: 'prot_pool')
%
% Output:
%   ecModel    - Model with updated Umin and Xi

    if nargin < 5 || isempty(prot_pool)
        prot_pool = 'prot_pool';
    end

    % Find basal branch reaction
    basal_rxn = [ex_rxn '_basal'];
    basalIdx = find(strcmp(ecModel.rxns, basal_rxn));

    if isempty(basalIdx)
        warning('Basal reaction %s not found. Apply sluice structure first.', basal_rxn);
        return;
    end

    % Update basal branch bounds (Umin)
    ecModel.lb(basalIdx) = -umin;
    ecModel.ub(basalIdx) = 1000;

    % Find extended branch reaction
    ext_rxn = [ex_rxn '_extended'];
    extIdx = find(strcmp(ecModel.rxns, ext_rxn));

    if ~isempty(extIdx)
        % Extended branch: uptake beyond Umin with Xi coefficient
        ecModel.lb(extIdx) = -1000;
        ecModel.ub(extIdx) = -umin;

        % Update Xi coefficient in S matrix
        protRow = find(strcmp(ecModel.mets, prot_pool));
        if ~isempty(protRow)
            ecModel.S(protRow, extIdx) = xi;
        end
    end

    fprintf('[setSluiceParams] Set %s: Umin=%.4f, Xi=%.4f\n', ex_rxn, umin, xi);
end
