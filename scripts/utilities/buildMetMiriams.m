function metMiriams = buildMetMiriams(metAnnotation)
% BUILDMETMIRIAMS Constructs MIRIAM annotation structures for metabolites.
%
% Input:
%   metAnnotation - A cell array of metabolite annotation structures.
%                   Each cell contains a structure with various annotation fields.
%
% Output:
%   metMiriams   - A cell array where each element is a scalar structure 
%                  with two fields:
%                  'name'  - a cell array of annotation field names (with underscores replaced by dots),
%                  'value' - a cell array of the corresponding annotation values.

    numMets = numel(metAnnotation);
    metMiriams = cell(numMets, 1);
    
    for i = 1:numMets
        % Initialize an empty structure with fields 'name' and 'value'
        miriamStruct.name = {};
        miriamStruct.value = {};
        
        % Process if the current metabolite annotation is a structure.
        if isstruct(metAnnotation{i})
            % Extract all annotation field names.
            annFields = fieldnames(metAnnotation{i});
            % Initialize new cell arrays for expanded field names and values.
            newAnnFields = {};
            newAnnValues = {};
            % Loop over each annotation field.
            for j = 1:length(annFields)
                currentField = annFields{j};
                currentValue = metAnnotation{i}.(currentField);
                if iscell(currentValue)
                    % If the value is a cell array, replicate the field name for each element.
                    newAnnFields = [newAnnFields; repmat({currentField}, numel(currentValue), 1)];
                    % Flatten the cell array and append each element.
                    newAnnValues = [newAnnValues; currentValue(:)];
                else
                    % Otherwise, add the field name and its value once.
                    newAnnFields = [newAnnFields; {currentField}];
                    newAnnValues = [newAnnValues; {currentValue}];
                end
            end
            % Replace underscores with dots in newAnnFields.
            newAnnFields = cellfun(@(s) strrep(s, '_', '.'), newAnnFields, 'UniformOutput', false);
            % Build the scalar structure with fields 'name' and 'value'.
            miriamStruct.name = newAnnFields;
            miriamStruct.value = newAnnValues;
        end
        metMiriams{i} = miriamStruct;
    end
end