function model = Parameterization(model, varargin)
% PARAMETERIZATION Updates model parameters via a single dialog box.
% This function prompts the user to update the following fields:
%   Biological context:
%       - Organism
%       - Taxonomic ID
%       - UniProt Identifier System
%       - UniProt IDs
%   Enzyme constraints:
%       - Total protein content [g/gDW] (Ptot)
%       - Mass fraction of enzymes (f)
%       - Enzyme saturation factor (sigma)
%
% If the user cancels, the original model is returned unchanged.

% If using dialog, canceling returns the original model.
% Save original model in case of cancel
orig_model = model;

% Check for name-value pair inputs
if nargin > 1
    % Validate even number of additional arguments
    if mod(numel(varargin), 2) ~= 0
        error('Additional arguments must be provided as name-value pairs.');
    end
    
    % Process each parameter
    for i = 1:2:numel(varargin)
        paramName = varargin{i};
        value = varargin{i+1};
        
        switch lower(paramName)
            case 'organism'
                model.information.organism = convertToString(value);
            case 'taxonomicid'
                model.information.taxonomicID = convertToString(value);
            case 'uniprot_type'
                model.information.uniprot_type = convertToString(value);
            case 'uniprot_id'
                model.information.uniprot_id = convertToString(value);
            case 'uniprot_geneidfield'
                model.information.uniprot_geneidfield = convertToString(value);
            case 'uniprot_reviewed'
                model.information.uniprot_reviewed = value;
            case 'ptot'
                model.enzymeConstraints.Ptot = convertToNumber(value);
            case 'f'
                model.enzymeConstraints.f = convertToNumber(value);
            case 'sigma'
                model.enzymeConstraints.sigma = convertToNumber(value);
            otherwise
                warning('Parameter "%s" is not recognized and will be ignored.', paramName);
        end
    end
else
    % Dialog-based input
    prompt = { ...
        'Organism:', ...
        'Taxonomic ID:', ...
        'UniProt Identifier System:', ...
        'UniProt ID:', ...
        'UniProt geneIDfield:', ...
        'UniProt reviewed:', ...
        'Total protein content [g/gDW] (Ptot):', ...
        'Mass fraction of enzymes (f):', ...
        'Enzyme saturation factor (sigma):'};
    
    dlg_title = 'Model Parameterization';
    num_lines = 1;
    
    defaultans = { ...
        num2str(model.information.organism), ...
        num2str(model.information.taxonomicID), ...
        num2str(model.information.uniprot_type), ...
        num2str(model.information.uniprot_id), ...
        num2str(model.information.uniprot_geneidfield), ...
        num2str(model.information.uniprot_reviewed), ...
        num2str(model.enzymeConstraints.Ptot), ...
        num2str(model.enzymeConstraints.f), ...
        num2str(model.enzymeConstraints.sigma)};
    
    answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
    
    if isempty(answer)
        model = orig_model;
        return;
    end
    
    % Update model with dialog inputs
    model.information.organism = answer{1};
    model.information.taxonomicID = answer{2};
    model.information.uniprot_type = answer{3};
    model.information.uniprot_id = answer{4};
    model.information.uniprot_geneidfield = answer{5};
    model.information.uniprot_reviewed = answer{6};
    model.enzymeConstraints.Ptot = str2double(answer{7});
    model.enzymeConstraints.f = str2double(answer{8});
    model.enzymeConstraints.sigma = str2double(answer{9});
end
end

% Helper functions for type conversion
function str = convertToString(value)
    if isnumeric(value)
        str = num2str(value);
    else
        str = char(value);
    end
end

function num = convertToNumber(value)
    if ischar(value)
        num = str2double(value);
    else
        num = double(value);
    end
end
