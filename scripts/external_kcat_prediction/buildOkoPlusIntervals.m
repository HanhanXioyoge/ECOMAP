function result = buildOkoPlusIntervals(model, varargin)
%BUILDOKOPLUSINTERVALS Public MATLAB API for cross-species OKO+ intervals.
%
%   result = buildOkoPlusIntervals(model)
%   result = buildOkoPlusIntervals(model, optionsStruct)
%   result = buildOkoPlusIntervals(model, 'Predictors', {'UniKP'}, ...)
%
%   The function returns a plain MATLAB struct (no web/JSON envelope).  Web
%   adapters call the same core and add their transport envelope separately.

    options = parseInputs(varargin{:});
    predictors = take(options, 'Predictors', {'UniKP'});
    parameters = take(options, 'Parameters', []);
    options = rmfieldIfPresent(options, {'Predictors','Parameters'});
    result = runOkoPipeline(model, predictors, options, parameters);
end

function options = parseInputs(varargin)
    options = struct();
    if isempty(varargin), return; end
    if isscalar(varargin) && isstruct(varargin{1})
        options = varargin{1};
        return;
    end
    if mod(numel(varargin), 2) ~= 0
        error('buildOkoPlusIntervals:InvalidOptions', ...
            'Options must be a struct or name-value pairs.');
    end
    for i = 1:2:numel(varargin)
        options.(char(varargin{i})) = varargin{i+1};
    end
end

function value = take(s, name, fallback)
    if isfield(s, name) && ~isempty(s.(name)), value = s.(name); else, value = fallback; end
end

function s = rmfieldIfPresent(s, names)
    for i = 1:numel(names)
        if isfield(s, names{i}), s = rmfield(s, names{i}); end
    end
end
