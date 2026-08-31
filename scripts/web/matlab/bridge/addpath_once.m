function addpath_once(dirPath)
%ADDPATH_ONCE Add a directory to MATLAB path the first time, no-op thereafter.
% Uses a persistent set keyed on the absolute, canonical path.
    persistent SEEN
    if isempty(SEEN), SEEN = {}; end
    canon = char(java.io.File(dirPath).getCanonicalPath());  % requires JVM
    for i = 1:numel(SEEN)
        if strcmp(SEEN{i}, canon), return; end
    end
    addpath(dirPath);
    SEEN{end+1} = canon;
end