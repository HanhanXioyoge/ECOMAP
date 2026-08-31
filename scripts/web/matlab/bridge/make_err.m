function payload = make_err(error_code, error_message)
%MAKE_ERR Build a failure envelope. error_code must be in the 14-code catalogue
% (validated against a hard-coded cell array; throw if unknown).
    if nargin < 2 || isempty(error_message)
        error_message = error_code;
    end
    % Validate against catalogue.
    catalogue = bridge_error_catalogue();
    if ~any(strcmp(error_code, catalogue))
        error('bridge:UnknownCode', ...
              'make_err: "%s" is not in the 14-code catalogue. Update CONTRACT.md first.', error_code);
    end
    payload = struct('ok', false, ...
                     'error_code', error_code, ...
                     'error_message', error_message, ...
                     'result', []);
end