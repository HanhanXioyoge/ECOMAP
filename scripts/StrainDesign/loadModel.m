function model = loadModel(model_path, model_type)
% scripts/StrainDesign/loadModel.m
% StrainDesign.loadModel  Thin wrapper around Reconstruction/loadModel that
% additionally marks the result with `information.ecModel = sniffModelType(model)`.
%
%   model = scripts.StrainDesign.loadModel(model_path)
%   model = scripts.StrainDesign.loadModel(model_path, model_type)
%
% Inputs:
%   model_path - char or string, path to the model file (.xml / .json / .yml).
%   model_type - char or string, optional. Defaults to 'Tradition' if empty/missing.
%
% Output:
%   model - struct, whatever Reconstruction/loadModel returns, augmented with
%           model.information.ecModel = 'ecGEM' | 'GEM'  (sniffModelType result;
%           'ecGEM' when the model is an enzyme-constrained GEM, e.g. has
%           enzyme-usage or proteomics pool reactions; 'GEM' for plain GEMs).
%
% Notes:
%   - Forward compatible with ECOMAP 2.0 strain-design tracks where loading a
%     model always implies the answer to "is this an ecModel?".
%   - Sniffing is delegated to scripts/StrainDesign/sniffModelType.m (migrated
%     from MDP; logic is heuristic and does not require any external toolbox).
%
% See also: scripts.Reconstruction.loadModel, scripts.StrainDesign.sniffModelType.

    if nargin < 2 || isempty(model_type)
        model_type = 'Tradition';
    end

    % Make Reconstruction reachable for the recursive call below.
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'Reconstruction'));

    % Delegate. Reconstruction/loadModel may need (filename, modeltype, modelDir, parameters).
    % Passing only 2 args is fine - the other two resolve from defaults inside it.
    % Unqualified recursion: MATLAB resolves to whichever loadModel is currently
    % ahead on the path - and we just put Reconstruction/ first.
    model = loadModel(char(model_path), char(model_type));

    % Annotate with the ecModel flag the strain-design pipeline consumes.
    if ~isfield(model, 'information') || ~isstruct(model.information)
        model.information = struct();
    end
    model.information.ecModel = sniffModelType(model);
end