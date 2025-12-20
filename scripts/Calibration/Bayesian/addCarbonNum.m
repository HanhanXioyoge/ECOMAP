function model = addCarbonNum(model,bioRxn)
% addCarbonNum
%   This function calculates the carbon number for each exchange reaction and adds it to the model.
%   It then computes the carbon content for the biomass reaction using the function `getcarbonnum_for_growth`.

% 1. Initialize carbon content for all reactions
model.excarbon = zeros(1,length(model.rxns));

% 2. Get the exchange reactions and their indices
[EXrxn, EXrxnIdx] = getExchangeRxns(model);

% 3. Get the carbon numbers for each exchange reaction
[CarbonNum,EXfors] = getcarbonnum(model,EXrxn);

% 4. Store the carbon numbers for the exchange reactions in the model
model.excarbon(EXrxnIdx) = CarbonNum;

% 5. Get the carbon content for the biomass reaction and update the model
model.excarbon(strcmp(model.rxns,bioRxn)) = getcarbonnum_for_growth(model, bioRxn);

end

function [CarbonNum,EXfors,changedmets] = getcarbonnum(model,exrxn)
% getcarbonnum
%   This function calculates the carbon number for each exchange reaction.
%   It assumes that the input exchange reactions contain only one metabolite in the equation.

% 1. Find the index of the exchange reactions in the model
[~,idx] = ismember(exrxn,model.rxns);

% 2. Get the metabolites associated with the exchange reactions
EXmets = model.S(:,idx);

% 3. Find the metabolites involved in each exchange reaction
EXmetsIdx = zeros(length(exrxn),1);
for k = 1:length(EXmets(1,:))
    EXmetsIdx(k) = find(EXmets(:,k)); 
end

% 4. Get the formulas of the metabolites
EXfors = model.metFormulas(EXmetsIdx);

% 5. Calculate the elemental composition for carbon atoms (C) in each metabolite
[Ematrix, elements] = getElementalComposition(EXfors,{'C'});
Ematrix = Ematrix(:,1);  % Get the carbon counts

% 6. Identify metabolites without carbon content (invalid entries)
changedmets = exrxn(isnan(Ematrix));

% 7. Replace missing carbon counts with 1
Ematrix(isnan(Ematrix)) = 1;

% 8. Return the carbon numbers for each exchange reaction
CarbonNum = Ematrix';

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

% 1. Get the formulas of the metabolites
formulas = model.metFormulas(metIndices);

% 2. Get the elemental composition for carbon atoms (C) in the metabolites
[Ematrix, ~] = getElementalComposition(formulas, {'C'});
C_counts = Ematrix(:, 1);  % Carbon count for each metabolite

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
