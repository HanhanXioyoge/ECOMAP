function model = loadYMLmodel(model, filename)
    % Load YML model from file
    try
        yml_model = readYAMLmodel(filename);
    catch
        error('Required function for loading YAML models is missing or failed.');
    end
    
    % Get all fieldnames from the yml_model
    fields = fieldnames(yml_model);
    
    % Loop over each field in the YAML model and copy it to model if that field already exists in model
    for i = 1:length(fields)
        field = fields{i};
        if isfield(model, field)
            model.(field) = yml_model.(field);
        end
    end
    
    % Check if the yml_model contains the 'ec' field for enzyme constraints
    if isfield(yml_model, 'ec')
        % Determining ecModel type
        if ~isfield(yml_model, 'ec') || ~isfield(yml_model.ec, 'geckoLight')
            error('The field yml_model.ec.geckoLight does not exist.');
        % Assign ecModeltype based on geckoLight value
        elseif yml_model.ec.geckoLight == true
            model.enzymeConstraints.ecModeltype = 'simple';
        elseif yml_model.ec.geckoLight == false
            model.enzymeConstraints.ecModeltype = 'complex';
        else
            error('yml_model.ec.geckoLight must be logical.');
        end

        % Define the allowed fields to be assigned from yml_model.ec
        allowed_fields = {'rxns', 'kcat', 'source', 'notes', 'eccodes', 'concs', 'genes', 'enzymes', 'mw', 'sequence', 'rxnEnzMat'};
        % Get the fieldnames from the yml_model.ec structure
        ec_fields = fieldnames(yml_model.ec);
        % Iterate over each field in yml_model.ec
        for i = 1:length(ec_fields)
            field = ec_fields{i};
            % If the current field is in the allowed list, assign it to model.enzymeConstraints
            if ismember(field, allowed_fields)
                model.enzymeConstraints.(field) = yml_model.ec.(field);
            end
        end
    end
end