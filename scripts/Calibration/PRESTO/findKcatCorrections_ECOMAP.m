function [kcatOrig, kcatCorr, protIDs, reportTable] = findKcatCorrections_ECOMAP(origModel, corrModel, enzMetPfx)
% FINDKCATCORRECTIONS_ECOMAP
% Compare kcat values between original and corrected ecModel (ECOMAP style).
% For each EC reaction whose kcat has increased, return:
%   - original kcat
%   - corrected kcat
%   - reaction ID
%   - a single cell containing all involved enzyme IDs joined by commas.
%
% Input
%   origModel   : ecModel before PRESTO correction (with enzymeConstraints)
%   corrModel   : ecModel after PRESTO correction
%   enzMetPfx   : (optional) prefix like 'prot_' if you想在 PROTEIN ID 中保留/去掉
%
% Output
%   kcatOrig    : vector of original kcat values (only those that increased)
%   kcatCorr    : vector of corrected kcat values
%   protIDs     : cell array, each cell is 'enz1,enz2,...'
%   reportTable : summary table

    if nargin < 3
        enzMetPfx = '';
    end

    EC0 = origModel.enzymeConstraints;
    EC1 = corrModel.enzymeConstraints;

    kcatOrigAll = EC0.kcat(:);
    kcatCorrAll = EC1.kcat(:);

    idxChanged = kcatCorrAll > kcatOrigAll & ...
                 isfinite(kcatOrigAll) & isfinite(kcatCorrAll);

    kcatOrig = kcatOrigAll(idxChanged);
    kcatCorr = kcatCorrAll(idxChanged);

    rxnIDsAll = EC0.rxns(:);
    rxnIDs    = rxnIDsAll(idxChanged);

    rxnEnzMat   = EC0.rxnEnzMat;
    allEnzNames = EC0.enzymes(:);

    changedRows = find(idxChanged);
    protIDs = cell(numel(changedRows), 1);

    for t = 1:numel(changedRows)
        r = changedRows(t);
        enzCols = find(rxnEnzMat(r, :) ~= 0); 

        enzList = allEnzNames(enzCols);

        if ~isempty(enzMetPfx)
            enzList = strcat(enzMetPfx, enzList);
        end

        protIDs{t} = strjoin(enzList, ',');
    end

    reportTable = unique( ...
        sortrows( ...
            table( ...
                protIDs, ... 
                rxnIDs, ...
                kcatOrig, ...
                kcatCorr, ...
                kcatCorr ./ kcatOrig, ...
                kcatCorr - kcatOrig, ...
                'VariableNames', { ...
                    'PROTEIN ID', ...
                    'REACTION ID', ...
                    'KCAT ORIG [s^-1]', ...
                    'KCAT UPDATED [s^-1]', ...
                    'FOLD-INCREASE', ...
                    'DELTA' ...
                }), ...
            'DELTA', 'descend'), ...
        'stable');

end
