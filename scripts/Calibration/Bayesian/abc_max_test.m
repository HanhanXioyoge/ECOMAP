function [rmse_final,exp,simulated,growthdata,UnconstrainedGrowth]=abc_max_test(ecModel,kcat_random_all,growthdata,UnconstrainedGrowth,proc,sample_generation,j,bioRxn,c_source,rxn2block,org_name)

nstep = sample_generation/proc;
rmse_final = zeros(1,nstep);
kcat_sample = kcat_random_all(:,(j-1)*nstep+1:j*nstep);

% get carbonnum for each exchange rxn to further calculation of error
if ~isfield(ecModel,'excarbon')
    ecModel = addCarbonNum(ecModel,bioRxn);
end

exp = [];
simulated = [];

for k = 1:nstep
    %disp(['nstep:',num2str(k),'/',num2str(nstep)])
    kcat_random  = kcat_sample(:,k);
    ecModel.enzymeConstraints.kcat = kcat_random;
    ecModel = UpdateSmatrix(ecModel);
    
    %% first search with substrate constraints
    objective = bioRxn;
    if ~isempty(growthdata)
        [rmse_1,exp_1,simulated_1] = rmsecal(ecModel,growthdata,true,objective,c_source,rxn2block,org_name);
    else
        rmse_1 = [];
        exp_1 = [];
        simulated_1 = [];
    end
    %% second search for maxmial growth rate without constraints
    if ~isempty(UnconstrainedGrowth)  % simulate the maximal growth rate
        [rmse_2,exp_2,simulated_2] = rmsecal(ecModel,UnconstrainedGrowth,false,objective,c_source,rxn2block,org_name);
    else
        rmse_2 = [];
        exp_2 = [];
        simulated_2 = [];
    end
    exp = [exp_1;exp_2];
    simulated = [simulated_1;simulated_2];
    rmse_final(1,k) = mean([rmse_1,rmse_2],'omitnan');
    
    %% only output simulated result for one generation
    if nstep ~= 1 || sample_generation ~= 1
        simulated = [];
        exp = [];
    end
end
end

function [rmse,exp_complete,simulated] = rmsecal(ecModel,data,constrain,objective,c_source,rxn2block,org_name)

    num_exps = height(data);
    simulated = zeros(num_exps, 9); 
    rmse_tmp = zeros(num_exps, 1);

    constraints = data.Properties.VariableNames;
    ex_mets          = constraints(4:10);
    [~,idx] = ismember(ex_mets,ecModel.rxns);
    exp_table = data(:, 4:10);
    exp = table2array(exp_table);  % ace eth gly pyr co2 o2 nh4
    exp_table = data(:, 2:10);
    exp_complete = table2array(exp_table);  % sub u ace eth gly pyr co2 o2 nh4

    for i = 1:num_exps
        model_tmp = ecModel;
        substrate_name = data{i, 1};
        substrate_uptake = data{i, 2};
        growth_rate = data{i, 3};

        model_tmp.c = double(strcmp(model_tmp.rxns, objective));

        if strcmp(data{i,11},'anaerobic') ||strcmp(data{i,11},'limited') 
            model_tmp = anaerobicModel(model_tmp,org_name);
        end

        if ~constrain
            nan_mask = isnan(exp(i,:));
            nonan_exp = exp(i, ~nan_mask)';
            idx_1 = idx(1:5);
            idx_2 = idx(6:7);
            model_tmp.lb(idx_1(nan_mask(1:5))) = 0;
            model_tmp.lb(idx_2(nan_mask(6:7))) = -1000;

            model_tmp.lb(idx(~nan_mask)) = nonan_exp;
            %No export of glucose
            model_tmp.lb(strcmp(model_tmp.rxns,c_source)) = 0; 

            model_tmp.lb(strcmp(model_tmp.rxns,objective)) = growth_rate;
            model_tmp.lb(strcmp(model_tmp.rxns,substrate_name)) = -1000;
            model_tmp.c = double(strcmp(model_tmp.rxns, substrate_name));

            sol_tmp = solveLP(model_tmp);
            if isempty(sol_tmp.f)
                exp_complete(i,1) = 999;

            % else
            %     exp_complete(i,1) = sol_tmp.f;
            end

            model_tmp.c = double(strcmp(model_tmp.rxns, objective));
            model_tmp.lb(strcmp(model_tmp.rxns,objective)) = 0;
            sol_tmp = solveLP(model_tmp);
        else
            %No export of glucose
            nan_mask = isnan(exp(i,:));
            nonan_exp = exp(i, ~nan_mask)';
            idx_1 = idx(1:5);
            idx_2 = idx(6:7);
            model_tmp.lb(idx_1(nan_mask(1:5))) = 0;
            model_tmp.lb(idx_2(nan_mask(6:7))) = -1000;

            model_tmp.lb(idx(~nan_mask)) = nonan_exp;
            %No export of glucose
            model_tmp.lb(strcmp(model_tmp.rxns,c_source)) = 0;
            if isnan(substrate_uptake)
                model_tmp.lb(strcmp(model_tmp.rxns,substrate_name)) = -1000;
            else
                model_tmp.lb(strcmp(model_tmp.rxns,substrate_name)) = substrate_uptake;
            end
            sol_tmp = solveLP(model_tmp);
        end
        

        if checkSolution(sol_tmp)
            sol(:,i) = sol_tmp.x;
            tmp = ~isnan(exp_complete(i,:));
            substrate_idx = find(strcmp(model_tmp.rxns, substrate_name));
            objective_idx = find(strcmp(model_tmp.rxns, objective));
            idx_complete = [substrate_idx, objective_idx, idx];
            
            excarbon = ecModel.excarbon(idx_complete);
            excarbon(excarbon == 0) = 1;
            exp_tmp = exp_complete(i,tmp).*excarbon(tmp);
            simulated_tmp = sol(idx_complete(tmp),i)'.*excarbon(tmp); % normalize the growth rate issue by factor 10
            
            exp_block = zeros(1,length(setdiff(rxn2block,model_tmp.rxns(idx_complete(2))))); % all zeros for blocked exchange mets exchange
            rxnblockidx = ismember(model_tmp.rxns,setdiff(rxn2block,model_tmp.rxns(idx_complete(2))));
            simulated_block = sol(rxnblockidx,i)'.* ecModel.excarbon(rxnblockidx); %
            exp_block = exp_block(simulated_block~=0);
            simulated_block = simulated_block(simulated_block~=0);
            if constrain
                rmse_tmp(i) = sqrt(immse([exp_tmp,exp_block], [simulated_tmp,simulated_block]));
            else
                if length(exp_tmp) >= 2
                    rmse_tmp(i) = sqrt(immse(exp_tmp(1:2), simulated_tmp(1:2)));
                else
                    rmse_tmp(i) = sqrt(immse(exp_tmp(1), simulated_tmp(1)));
                end

            end
            simulated(i,:) = sol(idx_complete,i)';
        else
            simulated(i,:) = NaN;
            rmse_tmp(i) = 999;
        end
    end
    rmse_tmp(isnan(rmse_tmp)) = []; %we just skip any case without solution, they are pretty rare, but exist
    rmse = sum(rmse_tmp)/height(data(:,1));
end
