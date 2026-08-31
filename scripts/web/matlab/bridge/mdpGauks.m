function result = mdpGauks(ec_model_id, gem_model_id, bio_rxn)
%MDPGAUKS Run GAUKS calibration.
%   result = mdpGauks(ec_model_id, gem_model_id, bio_rxn) returns the
%   standard bridge envelope (see CONTRACT.md) carrying:
%       .gauks_model_id  --  UUID of the calibrated ecModel
%       .summary_rows    --  row count of the GAUKS summary table
    addpath_once(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Calibration'));
    try
        ecModel = resolve_model_id(ec_model_id);
        gemModel = resolve_model_id(gem_model_id);
    catch err
        result = make_err('err_model_format', err.message);
        return;
    end
    if nargin < 3 || isempty(bio_rxn)
        result = make_err('err_no_biomass', 'bioRxn not provided');
        return;
    end
    bridge_log('mdpGauks', 'Running GAUKS');
    try
        [ecModel_out, ~, summaryTbl] = GAUKS(ecModel, gemModel, bio_rxn);
    catch err
        result = make_err('err_gurobi_license', err.message);
        return;
    end
    new_id = char(java.util.UUID.randomUUID.toString);
    register_model(new_id, ecModel_out);
    if isnumeric(summaryTbl)
        nrows = size(summaryTbl, 1);
    elseif isstruct(summaryTbl)
        nrows = numel(summaryTbl);
    elseif istable(summaryTbl)
        nrows = height(summaryTbl);
    else
        nrows = 0;
    end
    payload = struct( ...
        'gauks_model_id', new_id, ...
        'summary_rows', nrows);
    result = make_ok(payload);
end
