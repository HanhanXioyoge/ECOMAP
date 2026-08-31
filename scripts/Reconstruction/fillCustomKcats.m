function [model, rxnUpdated, notMatch] = fillCustomKcats(model, customKcats, rxnNameType, nameMap, parameters)
    % -------- Parameters --------
    % ECOMAP function to apply user-defined kcat values to the enzyme-constrained model.
    % This function updates kcat values for reactions based on the custom kcat data provided.
    % If customKcats is a file, it is read from the specified path. If the reaction names in the custom file 
    % differ from the model's reaction name format, it will replace them using the provided nameMap.
    %
    % Inputs:
    %   model        - ECOMAP model with enzyme constraints (contains enzymeConstraints structure)
    %   customKcats  - Either a structure containing custom kcat data or the path to a CSV file.
    %                  CSV columns should include: proteins, genes, gene_name, kcat, rxns, notes, stoicho.
    %   rxnNameType  - A string indicating the naming format of the reactions in the customKcats file.
    %   nameMap      - A table containing the mapping for reaction names based on ecModeltype.
    %                  Columns should represent different ecModeltype and rows should contain 
    %                  corresponding reaction names.
    %   parameters   - Optional parameter structure (if not provided, uses default parameters).
    %
    % Outputs:
    %   model        - Updated ECOMAP model with modified kcat values.
    %   rxnUpdated   - List of reaction IDs where kcat values have been updated.
    %   notMatch     - Table of reactions that could not be fully matched based on the custom kcat data.

    % Check and load parameters if not provided
    if nargin < 5 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end

    % Load custom kcat data if not provided
    if nargin < 2 || isempty(customKcats)
        customKcats = fullfile(parameters.reconstructionDir, 'kcatData', 'customKcats.csv');
    end
    
    % If the input is a structure, validate the essential fields
    if isstruct(customKcats)
        if ~all(isfield(customKcats, {'proteins', 'kcat', 'rxns'}))
            error('The customKcats structure does not have all essential fields.');
        end
    % If the input is a file, read it as a CSV
    elseif isfile(customKcats)
        data = readtable(customKcats, 'Delimiter', ',');  % Read CSV as table
        customKcats = struct();
        % Convert the table to a structure, with fields corresponding to the table columns
        customKcats.proteins = data.proteins;
        customKcats.genes = data.genes;
        customKcats.gene_name = data.gene_name;
        customKcats.kcat = data.kcat;
        customKcats.rxns = data.rxns;
        customKcats.notes = data.notes;
        customKcats.stoicho = data.stoicho;
    else
        error(['Cannot find file: ' customKcats]);
    end

    % Initialize reaction update flags and placeholders
    rxnToUpdate = false(length(model.enzymeConstraints.rxns), 1);
    rxnNotMatch = false(length(model.enzymeConstraints.rxns), 1);
    evaluatedRule = cell(length(model.enzymeConstraints.rxns), 1);
    enzInModel = cell(length(model.enzymeConstraints.rxns), 1);
    enzymeRxnNoSuffix = model.enzymeConstraints.rxns;

    % Process each custom kcat entry
    for i = 1:numel(customKcats.proteins)
        % Check if rxnNameType matches the model's enzymeConstraints.ecModeltype, replace names if needed
        if ~strcmp(rxnNameType, model.enzymeConstraints.ecModeltype)
            % Replace reaction names using nameMap based on rxnNameType
            if ismember(rxnNameType, nameMap.Properties.VariableNames)
                % Get the column corresponding to the rxnNameType and the model's ecModeltype
                nameColumn_custom = nameMap.(rxnNameType);
                nameColumn_model  = nameMap.(model.enzymeConstraints.ecModeltype);

                % Replace reactions in customKcats.rxns with the corresponding names from nameMap
                for j = 1:numel(customKcats.rxns)
                    rxnName = customKcats.rxns{j};
                    % Find the matching name in the custom column and replace it with the model's column
                    idx = find(strcmp(rxnName, nameColumn_custom));
                    if ~isempty(idx)
                        customKcats.rxns{j} = nameColumn_model{idx};  % Replace with the mapped name
                    end
                end
            else
                warning('The specified rxnNameType does not exist in the nameMap.');
            end
        end
        
        % If no specific proteins are provided, apply the kcat to the specified reactions
        if isempty(customKcats.proteins{i})
            rxns = strtrim(strsplit(customKcats.rxns{i}, ','));
            rxnIdxs = ismember(enzymeRxnNoSuffix, rxns);
            rxnToUpdate(rxnIdxs) = 1;
            model.enzymeConstraints.kcat(rxnIdxs) = customKcats.kcat(i);
        else
            % Process protein complexes (if multiple proteins are involved)
            prots = strtrim(strsplit(customKcats.proteins{i}, '+'));

            try
                enzIdx = cellfun(@(x) find(strcmpi(model.enzymeConstraints.enzymes, x)), prots);
            catch
                enzIdx = [];
                printOrange(['WARNING: Protein(s) ' customKcats.proteins{i} ' were not found in the model.']);
            end

            % If no reactions are specified, find all reactions that use the protein(s)
            if isempty(customKcats.rxns{i})
                temp_rxnIdxs = arrayfun(@(x) find(model.enzymeConstraints.rxnEnzMat(:, x)), enzIdx, 'UniformOutput', false);
            else
                rxns = strtrim(strsplit(customKcats.rxns{i}, ','));
                temp_rxnIdxs = arrayfun(@(x) find(strcmpi(enzymeRxnNoSuffix, x)), rxns, 'UniformOutput', false);
            end

            % If reaction indices are found, process them
            if ~isempty(temp_rxnIdxs)
                rxnIdxs = [];
                for j = 1:numel(temp_rxnIdxs)
                    rxnIdxs = [rxnIdxs; temp_rxnIdxs{j}];
                end

                % Remove duplicates and proceed with updating kcat values
                rxnIdxs = unique(rxnIdxs);

                for j = 1:numel(rxnIdxs)
                    % Get all enzymes involved in the reaction
                    allEnzInRxn = find(model.enzymeConstraints.rxnEnzMat(rxnIdxs(j), :));

                    % Calculate the matching percentage
                    C = intersect(enzIdx, allEnzInRxn);
                    if numel(C) == numel(enzIdx) && numel(C) == numel(allEnzInRxn)
                        match = 1;
                    else
                        if numel(enzIdx) < numel(allEnzInRxn)
                            match = numel(C) / numel(allEnzInRxn);
                        else
                            match = numel(C) / numel(enzIdx);
                        end
                    end

                    % Update the kcat if the match is perfect or above 50%
                    if match == 1
                        rxnToUpdate(rxnIdxs(j)) = 1;
                        model.enzymeConstraints.kcat(rxnIdxs(j)) = customKcats.kcat(i);

                        % Add a note indicating the custom kcat application
                        model.enzymeConstraints.source{rxnIdxs(j), 1} = 'custom';
                        if isfield(customKcats, 'notes')
                            if isempty(model.enzymeConstraints.notes)
                                model.enzymeConstraints.notes = cell(size(model.enzymeConstraints.kcat));
                            end
                            if isempty(model.enzymeConstraints.notes{rxnIdxs(j), 1}) && ~isempty(customKcats.notes{i})
                                model.enzymeConstraints.notes{rxnIdxs(j), 1} = customKcats.notes{i};
                            else
                                model.enzymeConstraints.notes{rxnIdxs(j), 1} = [model.enzymeConstraints.notes{rxnIdxs(j), 1} ', ' customKcats.notes{i}];
                            end
                        end
                    elseif match >= 0.5 && match < 1
                        rxnNotMatch(rxnIdxs(j)) = 1;
                        evaluatedRule{rxnIdxs(j), 1} = customKcats.proteins{i};
                        enzInModel{rxnIdxs(j), 1} = strjoin(model.enzymeConstraints.enzymes(allEnzInRxn), ' + ');
                    end
                end
            end
        end
    end

    % Collect the updated reactions
    rxnUpdated = model.enzymeConstraints.rxns(find(rxnToUpdate));

    % Remove from the unmatched reactions those that have been updated
    remove = and(rxnToUpdate, rxnNotMatch);
    rxnNotMatch(remove) = 0;
    evaluatedRule(remove) = '';
    enzInModel(remove) = '';

    % Generate a table for unmatched reactions
    idRxns = model.enzymeConstraints.rxns(find(rxnNotMatch));
    if strcmp(model.enzymeConstraints.ecModeltype, 'basic')
        idRxns = regexprep(idRxns, '^[0-9]{3}_', '');
    end
    fullIdx = cellfun(@(x) find(strcmpi(model.rxns, x)), idRxns);
    rxnsNames = model.rxnNames(fullIdx);
    evaluatedRule = evaluatedRule(~cellfun('isempty', evaluatedRule));
    enzInModel = enzInModel(~cellfun('isempty', enzInModel));
    rules = model.grRules(fullIdx);
    notMatch = table(idRxns, rxnsNames, evaluatedRule, enzInModel, rules, ...
        'VariableNames', {'rxns', 'name', 'custom enzymes', 'enzymes in model', 'rules'});
end
