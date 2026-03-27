function ecModel = setUmin(ecModel, ex_rxn, umin, prot_pool)
% setUmin
%   Sets the Umin value for a sluice-structured exchange reaction.
%
% Input:
%   ecModel     - Model with sluice structure (from applySluiceStructure)
%   ex_rxn     - Exchange reaction name (e.g., 'EX_glc__D_e')
%   umin       - Minimum uptake rate (mmol/gDW/h)
%   prot_pool   - Protein pool name (default: 'prot_pool')
%
% Output:
%   ecModel    - Model with updated Umin

    if nargin < 4 || isempty(prot_pool)
        prot_pool = 'prot_pool';
    end

    % Find basal branch reaction
    basal_rxn = [ex_rxn '_basal'];
    basalIdx = find(strcmp(ecModel.rxns, basal_rxn));

    if isempty(basalIdx)
        warning('Basal reaction %s not found. Apply sluice structure first.', basal_rxn);
        return;
    end

    % Update basal branch bounds
    ecModel.lb(basalIdx) = -umin;
    % ecModel.ub(basalIdx) = 1000;

    %{
    % Find extended branch reaction
    ext_rxn = [ex_rxn '_extended'];
    extIdx = find(strcmp(ecModel.rxns, ext_rxn));

    if ~isempty(extIdx)
        % Extended branch: uptake beyond Umin
        ecModel.lb(extIdx) = -1000;
        ecModel.ub(extIdx) = -umin;  % Negative ub means beyond umin
    end
    %}

    fprintf('[setUmin] Set %s umin = %.4f\n', ex_rxn, umin);
end
