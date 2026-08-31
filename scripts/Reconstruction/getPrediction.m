function kcatList_prediction = getPrediction(model, DLmodel, parameters)
% getPrediction
%   Read predicted kcat CSV from reconstructionDir for one of {'DLKcat','UniKP','CatPred'}
%   and return a standardized kcatList struct compatible with selectKcatValue.
%
% Inputs:
%   model      :  ECOMAP ecGEM
%   DLmodel    : 'DLKcat' | 'UniKP' | 'CatPred'
%   parameters : struct with field reconstructionDir (folder containing the CSV files)
%
% Output:
%   kcatList_prediction : struct with fields
%       source     : 'DLKcat' | 'UniKP' | 'CatPred'
%       rxns       : cellstr, reaction IDs (must be in model.enzymeConstraints.rxns)
%       genes      : cellstr, gene IDs
%       substrates : cellstr, substrate names (must be in model.metNames)
%       kcats      : double (/s)

    % ------------------------- Parameter defaults -------------------------
    if nargin < 3 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters)
            error('ParameterManager is not set.');
        end
    end
    if ~isfield(parameters,'reconstructionDir') || isempty(parameters.reconstructionDir)
        error('parameters.reconstructionDir is required and must point to the folder of CSVs.');
    end
    reconstructionDir = fullfile(parameters.reconstructionDir, 'kcatData');

    if nargin < 2 || isempty(DLmodel)
        error('You must specify DLmodel as one of {''DLKcat'',''UniKP'',''CatPred''}.');
    end
    if isstring(DLmodel), DLmodel = char(DLmodel); end
    validModels = {'DLKcat','UniKP','CatPred'};
    if ~ismember(DLmodel, validModels)
        error('DLmodel must be one of {''DLKcat'',''UniKP'',''CatPred''}.');
    end

    % ------------------------- Resolve file path --------------------------
    inFile = fullfile(reconstructionDir, [DLmodel, '.csv']);
    if ~exist(inFile,'file')
        error('Prediction file not found: %s', inFile);
    end

    % ------------------------- Read CSV (robustly) ------------------------
    T = readtable(inFile, 'VariableNamingRule', 'preserve');

    % ------------------------- Column mapping by model --------------------
    switch DLmodel
        case {'DLKcat','UniKP'}
            reqCols = {'ReactionName','GeneID','Substrate_name','predicted_kcat'};
            rxnCol  = 'ReactionName';
            geneCol = 'GeneID';
            subCol  = 'Substrate_name';
            kcatCol = 'predicted_kcat';
        case 'CatPred'
            reqCols = {'ReactionName','GeneID','Substrate','Prediction_(s^(-1))'};
            rxnCol  = 'ReactionName';
            geneCol = 'GeneID';
            subCol  = 'Substrate';
            kcatCol = 'Prediction_(s^(-1))';
    end
    missing = setdiff(reqCols, T.Properties.VariableNames);
    if ~isempty(missing)
        error('[%s] Missing required columns: %s', DLmodel, strjoin(missing, ', '));
    end

    % ------------------------- Extract & normalize ------------------------
    rxns = ensureCellstr(T.(rxnCol));
    genes = ensureCellstr(T.(geneCol));
    subs  = ensureCellstr(T.(subCol));
    kcats = toDoubleVec(T.(kcatCol));  % tolerant to 'NA' / strings

    % ------------------------- Drop invalid rows --------------------------
    valid = isfinite(kcats);
    rxns = rxns(valid);
    genes = genes(valid);
    subs  = subs(valid);
    kcats = kcats(valid);

    if isempty(kcats)
        error('%s file has no numeric kcat values after filtering.', DLmodel);
    end

    % ------------------------- Sanity checks vs model ---------------------
    if ~isfield(model,'enzymeConstraints') || ~isfield(model.enzymeConstraints,'rxns')
        error('model.enzymeConstraints.rxns is required.');
    end
    matchMets = ismember(subs, model.metNames);
    if ~all(matchMets)
        bad = unique(subs(~matchMets));
        warning('[%s] Some substrates not found in model.metNames. Showing unique missing entries:', DLmodel);
        disp(bad);
        error('[%s] The prediction CSV likely does not match this model (substrates mismatch).', DLmodel);
    end

    matchRxns = ismember(rxns, model.enzymeConstraints.rxns);
    if ~all(matchRxns)
        bad = unique(rxns(~matchRxns));
        warning('[%s] Some reactions not found in model.enzymeConstraints.rxns. Showing unique missing entries:', DLmodel);
        disp(bad);
        % error('[%s] The prediction CSV likely does not match this model (reactions mismatch).', DLmodel);
    end

    % ------------------------- Build kcatList struct ----------------------
    kcatList_prediction = struct();
    kcatList_prediction.source     = DLmodel;
    kcatList_prediction.rxns       = rxns;
    kcatList_prediction.genes      = genes;
    kcatList_prediction.substrates = subs;
    kcatList_prediction.kcats      = kcats;
end

% ------------------------- helpers -------------------------
function out = ensureCellstr(x)
    % Convert string/cellstr/char/numeric to cellstr safely
    if iscell(x)
        if all(cellfun(@ischar, x))
            out = x(:);
        else
            out = cellfun(@toChar, x, 'UniformOutput', false);
            out = out(:);
        end
    elseif isstring(x)
        out = cellstr(x);
    elseif ischar(x)
        out = cellstr(x);
    elseif isnumeric(x)
        out = cellstr(string(x));
    else
        out = cellstr(string(x));
    end
end

function y = toChar(v)
    if isstring(v), y = char(v); 
    elseif ischar(v), y = v;
    elseif isnumeric(v), y = num2str(v);
    else, y = char(string(v));
    end
end

function v = toDoubleVec(col)
    % Convert kcat column to double vector; 'NA'/empty => NaN
    if iscell(col)
        v = nan(numel(col),1);
        for i = 1:numel(col)
            c = col{i};
            if isstring(c), c = char(c); end
            if ischar(c)
                c = strtrim(c);
                if isempty(c) || strcmpi(c,'NA')
                    v(i) = NaN;
                else
                    v(i) = str2double(c);
                end
            elseif isnumeric(c)
                v(i) = c;
            else
                v(i) = NaN;
            end
        end
    elseif isstring(col)
        v = str2double(col);
    elseif isnumeric(col)
        v = col(:);
    else
        v = str2double(string(col));
    end
end
