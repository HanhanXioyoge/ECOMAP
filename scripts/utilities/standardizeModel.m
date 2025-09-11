function model = standardizeModel(model, fileNameWithoutExtension, modeltype)
    % Core model identification
    model.id = fileNameWithoutExtension;          % Unique model identifier (string)
    model.name = modeltype;                       % Human-readable model name (string)
end