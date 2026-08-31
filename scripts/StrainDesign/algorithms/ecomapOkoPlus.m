function out = ecomapOkoPlus(model, biomassRxn, targetRxn, varargin)
%ECOMAPOKOPLUS Friendly alias for the public one-call OKO+ MATLAB API.
    out = runOkoPlus(model, biomassRxn, targetRxn, varargin{:});
end
