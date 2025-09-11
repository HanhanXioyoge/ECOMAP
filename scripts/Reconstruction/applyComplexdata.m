function [model, foundComplex, proposedComplex] = applyComplexdata(model, complexInfo, parameters)
% applyComplexdata  Map ComplexPortal complexes to model reactions and update stoichiometry
%
% Inputs
%   model         : ecModel having model.enzymeConstraints fields:
%                   rxns (cellstr), enzymes (cellstr), rxnEnzMat (double), ecModeltype (char)
%   complexInfo   : (a) struct array from ComplexPortal or
%                   (b) path to a JSON file produced by getComplexdata / ComplexPortal cache
%                   Expected fields per complex: complexID, name, geneName, protID, stochiometry
%   parameters    : (optional) parameter struct; if omitted, uses ParameterManager.getParams()
%
% Outputs
%   model            : updated model with rxnEnzMat filled with complex stoichiometry when exact match found
%   foundComplex     : table with exact matches (rxn, complexID, name, genes, protID_model, protID_complex, stochiometry)
%   proposedComplex  : table with suggested matches (same columns + match percentage)
%
% Note 
% Refer to the GECKO code

    % -------- Parameters and inputs --------
    if nargin < 3 || isempty(parameters)
        parameters = ParameterManager.getParams();
        if isempty(parameters), error('ParameterManager is not set.'); end
    end

    if nargin < 2 || isempty(complexInfo)
        complexInfo = getComplexdata(parameters.taxonomicID);
    end

    % Accept either a JSON path or already-decoded struct
    if ischar(complexInfo) || isstring(complexInfo)
        jsonStr     = fileread(complexInfo);
        complexData = jsondecode(jsonStr);
    else
        complexData = complexInfo;
    end

    if isempty(model.enzymeConstraints)
        error('The model does not contain the enzymeConstraints structure.');
    end
    ecModeltype = model.enzymeConstraints.ecModeltype;

    % Reaction names: strip prefix if needed
    switch lower(ecModeltype)
        case {'integrated','isozyme'}
            rxnNames = model.enzymeConstraints.rxns;
        otherwise
            rxnNames = extractAfter(model.enzymeConstraints.rxns, 4);
    end

    nRxn      = numel(rxnNames);
    enzNames  = model.enzymeConstraints.enzymes(:);
    rxnEnzMat = model.enzymeConstraints.rxnEnzMat;   % (#rxns x #enzymes)

    % -------- Normalize ComplexPortal data (protID & stochiometry) --------
    nCplx = numel(complexData);
    protCell   = cell(nCplx,1);  % each entry: column cellstr of protein IDs
    stoichCell = cell(nCplx,1);  % each entry: numeric row vector of stoichiometry

    for i = 1:nCplx
        % protID -> column cellstr
        p = complexData(i).protID;
        if iscell(p)
            if any(cellfun(@iscell, p)), p = horzcat(p{:}); end
            p = cellfun(@char, p, 'UniformOutput', false);
        elseif isstring(p) || ischar(p)
            p = {char(p)};
        else
            p = {char(string(p))};
        end
        if isrow(p), p = p(:); end
        protCell{i} = p;

        % stochiometry -> numeric row vector; all-zero -> ones
        s = complexData(i).stochiometry;
        if iscell(s)
            if all(cellfun(@isnumeric, s))
                s = cell2mat(s(:))';
            else
                s = strjoin(cellfun(@char, s, 'UniformOutput', false), ',');
            end
        end
        if isstring(s), s = char(s); end
        if ischar(s)
            s2 = strrep(s, ';', ',');                   % unify delimiter
            s2 = strtrim(regexprep(s2, '^\[|\]$', '')); % strip outer brackets if any
            toks = regexp(s2, '[, \t]+', 'split');
            nums = str2double(toks);
            if ~isempty(nums) && all(~isnan(nums))
                s = nums(:)';                            % row vector
            else
                s = zeros(1, numel(p));                  % fallback (align length)
            end
        end
        if isnumeric(s) && all(s == 0), s(:) = 1; end    % all-zero -> ones
        if iscolumn(s), s = s'; end

        % Length align between stoichiometry and protID
        if numel(s) < numel(p)
            s = [s, ones(1, numel(p) - numel(s))];
        elseif numel(s) > numel(p)
            s = s(1:numel(p));
        end
        stoichCell{i} = s;
    end

    % -------- Build complex-by-protein sparse matrix --------
    % Stack all proteins, uniquify, and create sparse indices
    allProt = vertcat(protCell{:});                           % long column of all protIDs
    [complexProts, ~, allProtCol] = unique(allProt, 'stable');
    countPerCplx = cellfun(@numel, protCell);
    rows = repelem((1:nCplx)', countPerCplx);                 % row indices per subunit
    cols = allProtCol;                                        % column indices per subunit
    vals = cell2mat(cellfun(@(v) v(:), stoichCell, 'UniformOutput', false));  % values

    complexMatrix = sparse(rows, cols, vals, nCplx, numel(complexProts));     % (#complex x #proteins)

    % -------- Pre-map model enzymes to complexMatrix columns --------
    [inComplexProt, enz2col] = ismember(enzNames, complexProts);

    % -------- Match complexes to reactions --------
    foundRows    = {};   % to accumulate exact matches
    proposedRows = {};   % to accumulate suggested matches

    % Optional progress printing (lightweight)
    % fprintf('Assigning complexes to %d reactions...\n', nRxn);
    
    for i = 1:nRxn
        enzMask = rxnEnzMat(i,:) ~= 0;                % enzymes used by this reaction
        if ~any(enzMask), continue; end
    
        modelProts = enzNames(enzMask);
        colIdx     = enz2col(enzMask);
        validMask  = inComplexProt(enzMask) & (colIdx > 0);
        if ~any(validMask), continue; end
    
        colIdx_valid     = colIdx(validMask);
        modelProts_valid = modelProts(validMask);
    
        % Candidate complexes: any nonzero in selected columns
        cand = any(complexMatrix(:, colIdx_valid) > 0, 2);
        if ~any(cand), continue; end
        potComplex = find(cand);

        % Stats: matched subunits, total subunits (in complex), model subunits
        subMatch   = sum(complexMatrix(potComplex, colIdx_valid) > 0, 2);
        totalUnits = sum(complexMatrix(potComplex, :) > 0, 2);
        modUnits   = numel(modelProts);
    
        percMatch  = subMatch ./ modUnits;       % coverage of model subunits
        relSize    = totalUnits ./ modUnits;     % complex size relative to model
    
        % ---- Exact match: model subunits == complex subunits, fully covered ----
        exactMask = (percMatch == 1) & (relSize == 1);
        if any(exactMask)
            ii = potComplex(find(exactMask, 1, 'first'));     % if multiple, take the first

            % Write back stoichiometry for involved enzymes (aligned order)
            [~, modelProtsIdx] = ismember(modelProts_valid, enzNames);
            valsAssign = full(complexMatrix(ii, colIdx_valid));
            rxnEnzMat(i, modelProtsIdx) = valsAssign;

            % Record exact match
            foundRows(end+1, :) = { ...
                model.enzymeConstraints.rxns{i}, ...
                complexData(ii).complexID, ...
                complexData(ii).name, ...
                complexData(ii).geneName, ...
                modelProts, ...
                complexData(ii).protID, ...
                complexData(ii).stochiometry ...
            }; %#ok<AGROW>

        % ---- Suggestions: only for multi-protein grRules ----
        elseif modUnits > 1
            % Case 1: model fully covered but complex contains more subunits
            moreMask = (percMatch == 1) & (relSize > 1);
            if any(moreMask)
                [~, k] = min(relSize(moreMask));             % choose minimal extra size
                jj = potComplex(find(moreMask, k, 'first')); % index subset
                jj = jj(k);                                  % pick selected
                proposedRows(end+1, :) = { ...
                    model.enzymeConstraints.rxns{i}, ...
                    complexData(jj).complexID, ...
                    complexData(jj).name, ...
                    complexData(jj).geneName, ...
                    modelProts, ...
                    complexData(jj).protID, ...
                    complexData(jj).stochiometry, ...
                    relSize(find(moreMask, 1, 'first') + k - 1) * 100 ...
                }; %#ok<AGROW>
            end

            % Case 2: >=75% coverage and complex not larger (comparable or smaller)
            partMask = (percMatch >= 0.75) & (percMatch < 1) & (relSize <= 1);
            if any(partMask)
                [~, k] = max(percMatch(partMask));           % choose highest coverage
                idxList = find(partMask);
                jj = potComplex(idxList(k));
                proposedRows(end+1, :) = { ...
                    model.enzymeConstraints.rxns{i}, ...
                    complexData(jj).complexID, ...
                    complexData(jj).name, ...
                    complexData(jj).geneName, ...
                    modelProts, ...
                    complexData(jj).protID, ...
                    complexData(jj).stochiometry, ...
                    percMatch(idxList(k)) * 100 ...
                }; %#ok<AGROW>
            end
        end
    end

    % Assign updated matrix back to the model
    model.enzymeConstraints.rxnEnzMat = rxnEnzMat;

    % -------- Output tables --------
    colNames = {'rxn','complexID','name','genes','protID_model','protID_complex','stochiometry'};
    if isempty(foundRows)
        foundComplex = cell2table(cell(0, numel(colNames)), 'VariableNames', colNames);
    else
        foundComplex = cell2table(foundRows, 'VariableNames', colNames);
    end

    if isempty(proposedRows)
        proposedComplex = cell2table(cell(0, numel(colNames)+1), 'VariableNames', [colNames, 'match']);
    else
        proposedComplex = cell2table(proposedRows, 'VariableNames', [colNames, 'match']);
    end

    fprintf('A total of %d complexes have full match, and %d proposed.\n', height(foundComplex), height(proposedComplex));
end
