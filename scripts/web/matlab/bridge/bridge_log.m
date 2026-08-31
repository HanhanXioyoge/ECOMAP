function bridge_log(name, varargin)
%BRIDGE_LOG Standard "[mdpName] msg" stdout line.
    fprintf(1, '[mdp%s] %s\n', name, sprintf(varargin{:}));
end