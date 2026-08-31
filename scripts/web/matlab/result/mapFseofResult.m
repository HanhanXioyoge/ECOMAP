function tbl = mapFseofResult(result)
% mapFseofResult  Adapt algFseof output (.rows) to the unified candidates table.
    arguments
        result struct
    end
    n = numel(result.rows);
    reaction     = cell(n,1);
    name         = cell(n,1);
    direction    = zeros(n,1);
    intervention = cell(n,1);
    score        = zeros(n,1);
    scoreLabel   = repmat({'FSEOF slope'}, n, 1);
    genes        = cell(n,1);
    for i = 1:n
        row = result.rows(i);
        reaction{i}     = char(row.enzymeID);
        name{i}         = char(row.enzymeName);
        direction(i)    = row.direction;
        score(i)        = row.slope;
        genes{i}        = char(row.grRule);
        if row.slope >= 0, intervention{i} = 'OE'; else, intervention{i} = 'KD'; end
    end
    tbl = table(reaction, name, direction, intervention, score, scoreLabel, genes);
end
