function modelType = sniffModelType(model)
% sniffModelType  Classify a model as 'GEM' or 'ecGEM' (GECKO 3 heuristic).
%   Score-based: enzymes, MWs, pathways, prot_pool_exchange. Conservative:
%   defaults to 'GEM' unless enzyme-constrained evidence is strong.
    arguments
        model struct
    end
    score = 0;
    if isfield(model, 'enzymes'),            score = score + 1; end
    if isfield(model, 'MWs'),                score = score + 1; end
    if isfield(model, 'pathways'),           score = score + 1; end
    if isfield(model, 'prot_pool_exchange'), score = score + 2; end
    if isfield(model, 'ec'),                 score = score + 2; end   % GECKO 3 model.ec
    if score >= 4
        modelType = 'ecGEM';
    else
        modelType = 'GEM';
    end
end
