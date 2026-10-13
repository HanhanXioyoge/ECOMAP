function model = addCarbonNum(model,bioRxn)
% addCarbonNum
%   This function calculates the carbon number for each exchange reaction and adds it to the model.
%   It then computes the carbon content for the biomass reaction using the function `getcarbonnum_for_growth`.

% 1. Initialize carbon content for all reactions
model.excarbon = zeros(1,length(model.rxns));

% Include named exchange reactions even after a sluice transformation. Such
% reactions may contain a virtual metabolite and are no longer recognised by
% exchange-reaction finders that require a single metabolite.
exchangeIdx = find(startsWith(model.rxns, 'EX_'));
singleMetIdx = find(full(sum(spones(model.S), 1)) == 1);
exchangeIdx = unique([exchangeIdx(:); singleMetIdx(:)]);

for k = exchangeIdx(:)'
    metIdx = find(model.S(:, k));
    carbonCounts = NaN(numel(metIdx), 1);
    for j = 1:numel(metIdx)
        carbonCounts(j) = carbonFromFormula(model.metFormulas{metIdx(j)});
    end
    carbonCounts = carbonCounts(isfinite(carbonCounts));
    if ~isempty(carbonCounts)
        % Sluice branches contain the same substrate on both sides, plus
        % protein pool. Use the substrate's atom count, not their sum.
        model.excarbon(k) = max(carbonCounts);
    end
end

% 5. Get the carbon content for the biomass reaction and update the model
model.excarbon(strcmp(model.rxns,bioRxn)) = getcarbonnum_for_growth(model, bioRxn);

end

function nCarbon = carbonFromFormula(formula)
% Count only the carbon element; e.g. Ca and Co must not count as carbon.
    if isempty(formula)
        nCarbon = NaN;
        return;
    end
    parts = regexp(char(formula), '([A-Z][a-z]?)([0-9]*\.?[0-9]*)', 'tokens');
    nCarbon = 0;
    for i = 1:numel(parts)
        if strcmp(parts{i}{1}, 'C')
            count = str2double(parts{i}{2});
            if isnan(count), count = 1; end
            nCarbon = nCarbon + count;
        end
    end
end

function C_Biomass = getcarbonnum_for_growth(model, bioRxnSink)
% getcarbonnum_for_growth
%   This function calculates the net carbon content for the biomass reaction based on overall reaction carbon conservation.
%   The formula used is: Net_C = Sum(-1 * coefficient * carbon number), excluding the biomass metabolite itself.

% 1. Find the index of the sink reaction and the biomass metabolite
sinkRxnIdx = find(strcmp(model.rxns, bioRxnSink));
if isempty(sinkRxnIdx), error('Sink reaction not found.'); end

bioRxn_S_col_sink = model.S(:, sinkRxnIdx);

% Find the biomass metabolite index (usually has a coefficient of -1)
biomassMetIdx = find(bioRxn_S_col_sink < 0); 
if isempty(biomassMetIdx), error('Biomass metabolite not found.'); end

% 2. Find the production reaction (BOF)
% Look for columns where the coefficient in the S matrix is > 0
producerRxnIdx = find(model.S(biomassMetIdx(1), :) > 0);
if isempty(producerRxnIdx), error('BOF (production reaction) not found.'); end
mainBioRxnIdx = producerRxnIdx(1);

% 3. Get all metabolites involved in the BOF (both reactants and products)
bioRxn_S_col = model.S(:, mainBioRxnIdx);
all_met_indices = find(bioRxn_S_col ~= 0); % Get all metabolites involved in the reaction

% *** Key point: Exclude the biomass metabolite itself ***
all_met_indices(all_met_indices == biomassMetIdx(1)) = []; 

% 4. Get the corresponding original coefficients for the metabolites
coeffs = full(bioRxn_S_col(all_met_indices));

% --- Call the recursive function to calculate the carbon content ---
% Passing -coeffs to account for:
%   - Reactants (coeffs < 0) will become positive -> Add carbon to the total
%   - Products (coeffs > 0) will become negative -> Subtract lost carbon
C_Biomass = recursiveCarbonCalculation(model, all_met_indices, -coeffs, 0);

% 5. Check the result for carbon content
if C_Biomass <= 0 || isnan(C_Biomass)
    warning('C_Biomass calculation is abnormal, using fallback value 42.0');
    C_Biomass = 42.0; 
end

end

%% Recursive helper function
function total_C = recursiveCarbonCalculation(model, metIndices, current_net_coeffs, depth)
% recursiveCarbonCalculation
%   This recursive function calculates the total carbon content based on the provided metabolite indices, coefficients, and current depth.
%   It processes both metabolites with carbon formulas and pseudometabolites (which do not have formulas).

% Initialize total carbon count
total_C = 0;

% Limit recursion depth to avoid infinite loops (or excessive recursion)
if depth > 5, return; end 

% 1. Get carbon counts from the metabolite formulas
formulas = model.metFormulas(metIndices);
C_counts = NaN(numel(metIndices), 1);
for i = 1:numel(metIndices)
    C_counts(i) = carbonFromFormula(formulas{i});
end

% 3. Loop over each metabolite to calculate its contribution to carbon content
for i = 1:length(metIndices)
    met = metIndices(i);
    net_coeff = current_net_coeffs(i);
    
    if ~isnan(C_counts(i))
        % Case A: Metabolite has a chemical formula, directly accumulate carbon
        total_C = total_C + net_coeff * C_counts(i);
    elseif isempty(formulas{i}) || isnan(C_counts(i))
        % Case B: Pseudometabolite (no formula), dig deeper into its production reaction
        % Find the unique reaction producing the pseudometabolite (S_ij > 0)
        prodRxnIdx = find(model.S(met, :) > 0, 1); 
        
        if ~isempty(prodRxnIdx)
            sub_S_col = model.S(:, prodRxnIdx);
            % Find all components involved in the production reaction (both reactants and products)
            other_components = find(sub_S_col ~= 0);
            % Exclude the current pseudometabolite being broken down
            other_components(other_components == met) = [];
            
            % Get the output coefficient of the pseudometabolite in the production reaction (usually 1.0)
            product_out_coeff = model.S(met, prodRxnIdx);
            
            % Calculate the net coefficients for the next recursion level:
            % Subnet coefficients = current demand * (-1 * coefficients in the reaction / product output coefficient)
            sub_net_coeffs = net_coeff * (-sub_S_col(other_components) / product_out_coeff);
            
            % Recursively accumulate carbon from the next level
            total_C = total_C + recursiveCarbonCalculation(model, other_components, sub_net_coeffs, depth + 1);
        end
    end
end

end
