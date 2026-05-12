function ecModel = addEnzymeConstraintsToECMpy(ecModel)
% ADDENZYMECONSTRAINTSTOECMPY Add protein pool and enzyme constraints to ECMpy model

    if ~isfield(ecModel, 'enzymeConstraints')
        error('ecModel must have enzymeConstraints field');
    end

    EC = ecModel.enzymeConstraints;
    poolMetID = 'prot_pool';

    % --- Step 1: Add prot_pool metabolite directly ---
    if ~ismember(poolMetID, ecModel.mets)
        ecModel.mets{end+1, 1} = poolMetID;
        ecModel.metNames{end+1, 1} = poolMetID;
        ecModel.metComps{end+1, 1} = 'c';
        ecModel.metFormulas{end+1, 1} = '';
        ecModel.metCharges{end+1, 1} = NaN;
        ecModel.metMiriams{end+1, 1} = struct('name', {}, 'value', {});
        ecModel.b(end+1, 1) = 0;
        % Extend S matrix with new row
        ecModel.S(end+1, :) = 0;
    end

    % --- Step 2: Add prot_pool_exchange reaction using addRxns ---
    poolExID = 'prot_pool_exchange';
    if ~ismember(poolExID, ecModel.rxns)
        poolRxn.rxns = {poolExID};
        poolRxn.rxnNames = {poolExID};
        poolRxn.mets = {{poolMetID}};
        poolRxn.stoichCoeffs = {-1};
        poolRxn.lb = -1000;
        poolRxn.ub = 0;
        poolRxn.rev = 1;
        ecModel = addRxns(ecModel, poolRxn);
    end

    % --- Step 3: Find indices of enzyme-constrained reactions ---
    [found, rxnGlobalIdx] = ismember(EC.rxns, ecModel.rxns);

    % --- Step 4: Set coefficients in prot_pool row ---
    poolRowIdx = find(strcmp(ecModel.mets, poolMetID));
    validRxns = found & ~isnan(rxnGlobalIdx);

    for i = 1:numel(EC.rxns)
        if validRxns(i)
            rxnIdx = rxnGlobalIdx(i);
            kcat_MW_val = EC.kcat_MW(i);
            if isnan(kcat_MW_val) || kcat_MW_val == 0
                coeff = 0;
            else
                coeff = -1000 / kcat_MW_val;
            end
            ecModel.S(poolRowIdx, rxnIdx) = coeff;
        end
    end

    % --- Step 5: Set protein pool budget ---
    if isfield(EC, 'Ptot') && isfield(EC, 'f') && isfield(EC, 'sigma')
        Ptot = EC.Ptot;
        f = EC.f;
        sigma = EC.sigma;
        prot_pool_new = -(Ptot * f * sigma * 1000);
        ecModel.lb(strcmp(ecModel.rxns, poolExID)) = prot_pool_new;
    end
end