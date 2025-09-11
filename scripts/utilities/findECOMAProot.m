function [ecomapPath, prevDir] = findECOMAProot()
% findECOMAProot
%   Automatically locates the root directory of the ECOMAP framework by
%   traversing upward from the current function's location until it finds
%   the marker file 'LICENSE.md'.
%
% Outputs:
%   ecomapPath - Full path to the ECOMAP root directory
%   prevDir    - Directory where the function was called from (for return)

markerFile = 'LICENSE.md';  % Marker file used to identify ECOMAP root
ST = dbstack('-completenames');
prevDir = pwd();

% Get the path of this file (findECOMAProot.m)
thisFilePath = ST(strcmp({ST.name}, 'findECOMAProot')).file;
ecomapPath = fileparts(thisFilePath);  % Start from folder where this file is
rootFound = false;
while ~rootFound
    % Check if LICENSE.md exists in this directory
    if exist(fullfile(ecomapPath, markerFile), 'file') == 2
        rootFound = true;
    else
        parentPath = fileparts(ecomapPath);
        if strcmp(parentPath, ecomapPath)
            error(['Cannot find ECOMAP root. Marker file "' markerFile '" not found.']);
        end
        ecomapPath = parentPath;  % Move up one level
    end
end
end