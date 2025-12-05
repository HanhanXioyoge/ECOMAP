function [condNames, E, expVal, nutrExch, P] = getconditions(model, addM9medium, unit, parameters)

    % -------- Parameters --------
    if nargin < 4 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end
    dataDir = parameters.dataDir;
    fprintf('Reading experimental data from file...\n')

    % enzyme abundances
    protAbundanceData = readtable(fullfile(dataDir,'abs_proteomics.tsv'),...
        'FileType', 'text', 'ReadRowNames', true);
    condNames = protAbundanceData.Properties.VariableNames;

    if ~isempty(model)
        enzRxnIdx = find(contains(model.rxns, 'prot_') & ~ismember(model.rxns, {'prot_pool_exchange'}));
        enzMetIdx = find(any(model.S(:,enzRxnIdx),2) & ~ismember(model.mets, {'prot_pool'}));
    
        new_enzymes = cellfun(@(x)regexp(x,'[A-Z0-9]{6}','match'),model.metNames(enzMetIdx));
        if ~isempty(setdiff(new_enzymes, model.enzymeConstraints.enzymes))||~isempty(setdiff(model.enzymeConstraints.enzymes, new_enzymes))
            error('number of enzyme pseudometabolites differs from model.enzyme field. Unable to match proteomics data')
        elseif ~isequal(new_enzymes, model.enzymeConstraints.enzymes)
            warning(['Model.enzyme field does not have same ordering as enzyme pseudometabolites in stochiometric matrix\n' ...
                'Ordering E matrix according to pseudometabolites'])
            model.enzymeConstraints.enzymes = new_enzymes;
        end
        % get mws
        [~, loc]= ismember(new_enzymes, model.enzymeConstraints.enzymes);
        MWs  = model.enzymeConstraints.mw(loc);
        % match E to gene IDs in the model
        protIDs = protAbundanceData.Properties.RowNames;
        % map protein abundance data to proteins in the model
        idxModel2Data = cell2mat(cellfun(@(x)find(ismember(protIDs,x)),...
            model.enzymeConstraints.enzymes,'un',0));
        % set up enzyme abundance matrix
        E = nan(numel(model.enzymeConstraints.enzymes),size(protAbundanceData,2));
        E(ismember(model.enzymeConstraints.enzymes,protIDs),:) = table2array(protAbundanceData(idxModel2Data,:));
        E(isnan(E)) = 0;
        if strcmp(unit, 'mmol/gDW')
            E = E .*MWs;
        end
    else 
        E=[];
    end

    clear protAbundanceData
    
    % read experimental growth rates
    growthRateTab = readtable(fullfile(dataDir, 'growth_rates.tsv'), 'FileType', 'text', 'ReadRowNames', true);
    % if necessary, re-order according to condNames from protein abundance table
    matchIdx = cellfun(@(x)find(ismember(condNames,x)),growthRateTab.Properties.RowNames);
    expVal = table2array(growthRateTab(matchIdx,:));

    pTab = readtable(fullfile(dataDir, 'total_protein.tsv'), 'FileType','text','ReadRowNames',true);
    % if necessary, re-order according to condNames from protein abundance table
    matchIdx = cellfun(@(x)find(ismember(condNames,x)),pTab.Properties.RowNames);
    P = table2array(pTab(matchIdx,:));
    if ~all(P>0)
        warning('Total protein not available for all conditions')
        fprintf('==> Using maximum protein content for these conditions\n\n')
        P(isnan(P)) = max(P);
    end

    % nutrient uptake data
    nutrExch = readtable(fullfile(dataDir, 'csource.tsv'), 'FileType','text','ReadRowNames',true);
    % if necessary, re-order according to condNames from protein abundance table
    matchIdx = cellfun(@(x)find(ismember(condNames,x)),nutrExch.Properties.VariableNames);
    nutrExch = nutrExch(:,matchIdx);
    
    % add M9 medium
    % from doi: 10.1016/j.jbiotec.2009.10.007 (medium reference in Valgepea et
    % al.)
    % nutrients from MetaCyc M9 definition (w/o glycerol) and missing elements found by
    % biomassPrecursorCheck
    
    % medium used in Schmidt et al. 2015 further contains thiamine solution
    if addM9medium
        M9ExchRxns = {'EX_na1_e','EX_pi_e','EX_cl_e',...
            'EX_k_e','EX_nh4_e','EX_mg2_e','EX_so4_e','EX_mobd_e',...
            'EX_mn2_e','EX_ni2_e','EX_zn2_e','EX_cu2_e','EX_ca2_e',...
            'EX_fe2_e','EX_fe3_e','EX_cd2_e','EX_cobalt2_e','EX_h2o_e',...
            'EX_o2_e'};
        
        M9 = array2table(repelem(-1000,numel(M9ExchRxns),numel(condNames)),...
            'RowNames',M9ExchRxns,'VariableNames',nutrExch.Properties.VariableNames);
        
        nutrExch = [nutrExch;M9];
    end
end

