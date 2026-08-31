function intervals = buildOkoIntervals(predictions, varargin)
% buildOkoIntervals  Aggregate homolog predictions into OKO+ kcat ranges.
%
% Input rows represent predictions for homologous enzymes in different
% organisms and require reaction, enzyme, organism and kcat columns. Use:
%   T = buildOkoIntervals(predictions, 'Predictor','DLKcat')
%   T = buildOkoIntervals({dlk, unikp, catpred}, 'Consensus','intersection')
%
% Min/max is the paper-compatible legacy aggregation. Quantile aggregation
% is available for sensitivity analysis. Predictor-specific intervals are
% retained independently unless Consensus='intersection' is requested.

    ip=inputParser;
    addParameter(ip,'Predictor','');
    addParameter(ip,'Aggregation','minmax');
    addParameter(ip,'Quantiles',[0.05 0.95]);
    addParameter(ip,'MinimumHomologs',2);
    addParameter(ip,'Consensus','none');
    parse(ip,varargin{:}); o=ip.Results;
    if ~iscell(predictions), predictions={predictions}; end
    built=cell(numel(predictions),1);
    for i=1:numel(predictions)
        T=readPredictionTable(predictions{i});
        built{i}=aggregateOne(T,o);
    end
    intervals=vertcat(built{:});
    if strcmpi(o.Consensus,'intersection')
        intervals=intersectPredictors(intervals);
    elseif ~strcmpi(o.Consensus,'none')
        error('buildOkoIntervals:Consensus','Consensus must be none or intersection.');
    end
end

function out=aggregateOne(T,o)
    vars=lower(regexprep(string(T.Properties.VariableNames),'[^a-z0-9]',''));
    rx=findColumn(vars,["rxn","reaction","reactionid","rxnid","rxnname"]);
    enz=findColumn(vars,["uniprot","enzyme","enzymeid","protein","proteinid"]);
    org=findColumn(vars,["organism","species","taxon","taxid"]);
    kc=findColumn(vars,["kcat","predictedkcat","prediction","value"]);
    predictor=find(ismember(vars,["predictor","model","dlmodel"]),1);
    R=string(T{:,rx}); E=string(T{:,enz}); O=string(T{:,org}); K=toDouble(T{:,kc});
    valid=strlength(R)>0 & strlength(E)>0 & strlength(O)>0 & isfinite(K) & K>0;
    R=R(valid); E=E(valid); O=O(valid); K=K(valid);
    if isempty(o.Predictor)
        if isempty(predictor), predictorName="unspecified";
        else
            names=unique(string(T{valid,predictor}));
            if numel(names)~=1, error('buildOkoIntervals:MixedPredictors','Filter mixed predictors with Predictor.'); end
            predictorName=names;
        end
    else
        predictorName=string(o.Predictor);
        if ~isempty(predictor)
            keep=strcmpi(string(T{valid,predictor}),predictorName);
            R=R(keep);E=E(keep);O=O(keep);K=K(keep);
        end
    end
    [G,rxn,enzyme]=findgroups(R,E);
    rows=cell(max(G),1);
    for g=1:max(G)
        sel=G==g; [~,uniqueOrg]=unique(O(sel),'stable'); values=K(sel); values=values(uniqueOrg);
        if numel(values)<o.MinimumHomologs, continue; end
        if strcmpi(o.Aggregation,'minmax')
            bounds=[min(values),max(values)];
        elseif strcmpi(o.Aggregation,'quantile')
            bounds=quantile(values,o.Quantiles);
        else
            error('buildOkoIntervals:Aggregation','Aggregation must be minmax or quantile.');
        end
        rows{g}=table(rxn(g),enzyme(g),bounds(1),bounds(2),mean(values),median(values), ...
            numel(values),predictorName,string(o.Aggregation), ...
            'VariableNames',{'rxn','uniprot','min','max','mean','median','nHomologs','predictor','aggregation'});
    end
    out=vertcat(rows{~cellfun('isempty',rows)});
    if isempty(out)
        out=table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
            zeros(0,1),strings(0,1),strings(0,1),'VariableNames', ...
            {'rxn','uniprot','min','max','mean','median','nHomologs','predictor','aggregation'});
    end
end

function out=intersectPredictors(T)
    [G,rxn,enzyme]=findgroups(T.rxn,T.uniprot); rows=cell(max(G),1);
    for g=1:max(G)
        block=T(G==g,:); lo=max(block.min); hi=min(block.max);
        if height(block)<2 || lo>hi, continue; end
        rows{g}=table(rxn(g),enzyme(g),lo,hi,mean(block.mean),median(block.median), ...
            min(block.nHomologs),"consensus", "intersection", ...
            'VariableNames',T.Properties.VariableNames);
    end
    out=vertcat(rows{~cellfun('isempty',rows)});
end

function T=readPredictionTable(value)
    if istable(value), T=value;
    elseif ischar(value)||(isstring(value)&&isscalar(value)), T=readtable(char(value),'VariableNamingRule','preserve');
    else, error('buildOkoIntervals:Input','Each prediction input must be a table or delimited-file path.'); end
end
function idx=findColumn(vars,aliases)
    idx=find(ismember(vars,aliases),1);
    if isempty(idx), error('buildOkoIntervals:Columns','Missing required column: %s.',strjoin(aliases,', ')); end
end
function x=toDouble(x)
    if isnumeric(x),x=double(x);else,x=str2double(string(x));end
    x=x(:);
end
