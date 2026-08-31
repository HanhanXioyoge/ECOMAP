function out = runOkoPlus(model, biomassRxn, targetRxn, varargin)
%RUNOKOPLUS Public one-call MATLAB API for interval generation and solving.
%
%   out = runOkoPlus(model, biomassRxn, targetRxn, ...
%       'Predictors', {'UniKP'}, 'Profile', 'integrated')
%
%   Supply 'Intervals' as a table/CSV to skip interval generation.  When it
%   is omitted, buildOkoPlusIntervals is called and the first requested
%   predictor interval artifact is used.

    options = parseInputs(varargin{:});
    intervals = take(options, 'Intervals', []);
    profile = take(options, 'Profile', 'auto');
    buildResult = struct();

    if isempty(intervals)
        buildOptions = rmfieldIfPresent(options, solverOnlyFields());
        buildResult = buildOkoPlusIntervals(model, buildOptions);
        predictors = take(options, 'Predictors', {'UniKP'});
        if ~iscell(predictors), predictors = cellstr(predictors); end
        intervals = firstAvailableArtifact(buildResult.predictor_csv_paths, predictors);
    end

    solverOptions = take(options, 'SolverOptions', struct());
    solverOptions.Profile = profile;
    solution = algOkoPlus(model, biomassRxn, targetRxn, intervals, solverOptions);
    out = struct('solution', solution, 'intervalBuild', buildResult, ...
                 'intervals', intervals);
end

function path = firstAvailableArtifact(paths, predictors)
    for i = 1:numel(predictors)
        name = char(predictors{i});
        if isfield(paths, name) && ~isempty(paths.(name))
            path = paths.(name);
            return;
        end
    end
    error('runOkoPlus:NoIntervals', ...
        'No predictor produced a usable OKO+ interval artifact.');
end

function names = solverOnlyFields()
    names = {'Intervals','Profile','SolverOptions'};
end

function options = parseInputs(varargin)
    options = struct();
    if isempty(varargin), return; end
    if isscalar(varargin) && isstruct(varargin{1}), options = varargin{1}; return; end
    if mod(numel(varargin),2) ~= 0
        error('runOkoPlus:InvalidOptions','Options must be a struct or name-value pairs.');
    end
    for i=1:2:numel(varargin), options.(char(varargin{i}))=varargin{i+1}; end
end

function value = take(s,name,fallback)
    if isfield(s,name) && ~isempty(s.(name)), value=s.(name); else, value=fallback; end
end

function s = rmfieldIfPresent(s,names)
    for i=1:numel(names), if isfield(s,names{i}), s=rmfield(s,names{i}); end, end
end
