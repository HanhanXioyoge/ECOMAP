function growthRates = knockout_heatmap(ecModel, substrates, uptaken_rate, knockouts)
    % 输入:
    % ecModel: 酶约束代谢模型 (COBRA model structure)
    % substrates: 底物列表 (cell array of 24 substrates, e.g., {'glc__D_e', 'acetate_e', ...})
    % knockouts: 敲除列表 (cell array of genes to knock out, e.g., {'b0001', 'b0002', ...})
    
    % 初始化结果矩阵
    ecModel = setParam(ecModel, 'lb', 'EX_glc__D_e', 0);
    numKnockouts = length(knockouts);
    numSubstrates = length(substrates);
    growthRates = zeros(numSubstrates, numKnockouts);
    
    % 循环遍历每个敲除基因
    for i = 1:numKnockouts
        % 克隆原始模型
        knockoutModel = ecModel;
        
        % 敲除基因
        knockoutModel = knockOutGene(knockoutModel, knockouts{i});
        
        for j = 1:numSubstrates
            substrateExchangeRxn = substrates{j};
            mid_model = setParam(knockoutModel, 'lb', substrateExchangeRxn, uptaken_rate(j));
            % 优化模型并计算生物量
            sol = solveLP(mid_model);
            if isempty(sol.f)
                growthRates(j, i) = 0; 
            else
                growthRates(j, i) = sol.f;
            end
        end
    end
    
    % 绘制热图
    figure;
    heatmap(knockouts,substrates, growthRates, 'Colormap', jet, 'ColorbarVisible', 'on');
    title('Single Gene Knockouts Growth Rates on Different Substrates');
    xlabel('Knockouts');
    ylabel('Substrates');
end
