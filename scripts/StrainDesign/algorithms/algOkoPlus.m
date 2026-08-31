function result = algOkoPlus(model, biomassRxn, targetRxn, intervals, varargin)
% algOkoPlus  Run OKO+ with externally supplied cross-species kcat ranges.
%
%   intervals may be a table or a CSV path. Required semantic columns are
%   reaction, enzyme/UniProt, minimum kcat and maximum kcat (s^-1).
%   Predictor results must be supplied separately; this function does not
%   combine DLKcat, UniKP and CatPred point predictions into a false range.

    if nargin < 4 || isempty(intervals)
        error('algOkoPlus:MissingIntervals', ...
            'OKO+ requires a cross-species kcat interval table or CSV file.');
    end
    options = parseOptions(varargin{:});
    prepared = okoPrepareModel(model, biomassRxn, targetRxn, intervals, options);
    result = okoSolve(prepared, 'oko-plus', options);
end

function options = parseOptions(varargin)
    options = struct('Profile', 'auto');
    if isempty(varargin), return; end
    if isscalar(varargin) && isstruct(varargin{1})
        names = fieldnames(varargin{1});
        for i = 1:numel(names), options.(names{i}) = varargin{1}.(names{i}); end
        return
    end
    if mod(numel(varargin), 2) ~= 0
        error('algOkoPlus:InvalidOptions', 'Options must be a struct or name-value pairs.');
    end
    for i = 1:2:numel(varargin)
        options.(char(varargin{i})) = varargin{i+1};
    end
end
