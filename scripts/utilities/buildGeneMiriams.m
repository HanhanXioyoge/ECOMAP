function geneMiriams = buildGeneMiriams(genes)
% BUILDGENEMIRIAMS Constructs MIRIAM annotation structures for genes.
%
% Input:
%   genes - A cell array of gene structures. Each gene may have an 'annotation' field.
%
% Output:
%   geneMiriams - A cell array where each element is a scalar structure with two fields:
%                 'name' and 'value'. The 'name' field is a cell array of annotation field
%                 names, and the 'value' field is a cell array of the corresponding values.

    % Determine the number of genes
    numGenes = numel(genes);
    geneMiriams = cell(numGenes, 1);
    
    for i = 1:numGenes
        if isfield(genes{i}, 'annotation')
            % Extract annotation field names and their corresponding values
            annFields = fieldnames(genes{i}.annotation);
            annValues = struct2cell(genes{i}.annotation);
            % Build a scalar structure with fields 'name' and 'value'
            miriamStruct.name = annFields;
            miriamStruct.value = annValues;
        else
            % If there is no annotation field, assign empty cell arrays
            miriamStruct.name = {};
            miriamStruct.value = {};
        end
        geneMiriams{i} = miriamStruct;
    end
end