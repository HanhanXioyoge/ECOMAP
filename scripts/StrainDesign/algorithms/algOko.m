function result = algOko(model, biomassRxn, targetRxn, varargin)
% algOko  Minimize the number of kcat changes required for overproduction.
%
%   result = algOko(model, biomassRxn, targetRxn)
%   result = algOko(..., optionsStruct)
%   result = algOko(..., 'Name', value, ...)
%
% OKO follows Razaghi-Moghadam et al. (2024).  Profile='legacy' preserves
% the numerical conventions of the released scripts. Profile='integrated'
% reads ECOMAP enzymeConstraints and accounts for MW and subunit factors.

    options = parseOptions(varargin{:});
    prepared = okoPrepareModel(model, biomassRxn, targetRxn, [], options);
    result = okoSolve(prepared, 'oko', options);
end

function options = parseOptions(varargin)
    options = struct('Profile', 'auto');
    if isempty(varargin), return; end
    if isscalar(varargin) && isstruct(varargin{1})
        options = mergeStruct(options, varargin{1});
        return
    end
    if mod(numel(varargin), 2) ~= 0
        error('algOko:InvalidOptions', 'Options must be a struct or name-value pairs.');
    end
    for i = 1:2:numel(varargin)
        options.(char(varargin{i})) = varargin{i+1};
    end
end

function out = mergeStruct(out, extra)
    names = fieldnames(extra);
    for i = 1:numel(names), out.(names{i}) = extra.(names{i}); end
end
