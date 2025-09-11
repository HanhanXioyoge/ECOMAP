function model = loadXMLmodel(model, filename)
    % Load XML model
    try
        xml_model = importModel(filename);
    catch
        error('required for loading SBML models');
    end
    
    % Correspondence using names
    fields = fieldnames(xml_model);
    
    for i = 1:length(fields)
        field = fields{i};
        if isfield(model, field)
            model.(field) = xml_model.(field);
        end
    end
end