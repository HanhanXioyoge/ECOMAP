function rxnMiriams = buildRxnMiriams(rxns)
% BUILDRXNMIRIAMS Constructs MIRIAM annotation structures for reactions.
%
% Input:
%   rxns - A cell array of reaction structures. Each reaction may have an
%          'annotation' field.
%
% Output:
%   rxnMiriams - A cell array where each element is a scalar structure with two fields:
%                'name' and 'value'. For annotation fields that contain multiple elements,
%                the field name appears repeatedly (one for each element) and the values are
%                flattened into the 'value' cell array.

    numRxns = numel(rxns);
    rxnMiriams = cell(numRxns, 1);
    
    for i = 1:numRxns
        if isfield(rxns{i}, 'annotation')
            annFields = fieldnames(rxns{i}.annotation);   
            % Initialize new cell arrays for the expanded field names and values.
            newAnnFields = {};
            newAnnValues = {};
            % Loop over each annotation field.
            for j = 1:length(annFields)
                currentField = annFields{j};
                currentValue = rxns{i}.annotation.(currentField);
                if iscell(currentValue)
                    % If the value is a cell array, replicate the field name for each element.
                    newAnnFields = [newAnnFields; repmat({currentField}, numel(currentValue), 1)];
                    % Flatten the cell array and append each element.
                    newAnnValues = [newAnnValues; currentValue(:)];
                else
                    % Otherwise, just add the field name and its value once.
                    newAnnFields = [newAnnFields; {currentField}];
                    newAnnValues = [newAnnValues; {currentValue}];
                end
            end
            newAnnFields = cellfun(@(s) strrep(s, '_', '.'), newAnnFields, 'UniformOutput', false);
            % Build a scalar structure with fields 'name' and 'value'.
            miriamStruct.name = newAnnFields;
            miriamStruct.value = newAnnValues;
        else
            % If there is no annotation field, assign empty cell arrays.
            miriamStruct.name = {};
            miriamStruct.value = {};
        end
        rxnMiriams{i} = miriamStruct;
    end
end