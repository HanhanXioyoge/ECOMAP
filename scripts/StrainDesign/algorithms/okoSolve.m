function result = okoSolve(prepared, mode, options)
% okoSolve  Shared LP/MILP implementation for OKO and OKO+.
% Internal API: callers should normally use algOko or algOkoPlus.

    if nargin < 3, options = struct(); end
    mode = lower(char(mode));
    isPlus = strcmp(mode, 'oko-plus');
    if ~ismember(mode, {'oko','oko-plus'})
        error('okoSolve:UnknownMode', 'Mode must be oko or oko-plus.');
    end
    if isempty(which('gurobi'))
        error('okoSolve:GurobiRequired', ...
            'OKO requires the Gurobi MATLAB interface on the MATLAB path.');
    end

    o = defaults(options, prepared.profile);
    model = prepared.model;
    model.lb = double(model.lb(:)) * o.FluxScale;
    model.ub = double(model.ub(:)) * o.FluxScale;
    b = double(model.b(:)) * o.FluxScale;
    Aeq = sparse(double(model.S));
    if o.NormalizeRows
        rowMax = full(max(abs(Aeq), [], 2));
        nz = rowMax > 0 & isfinite(rowMax);
        Aeq(nz,:) = spdiags(1 ./ rowMax(nz), 0, nnz(nz), nnz(nz)) * Aeq(nz,:);
        b(nz) = b(nz) ./ rowMax(nz);
    end
    n = size(Aeq,2);
    bio = prepared.biomassIdx; target = prepared.targetIdx;
    params = struct('FeasibilityTol', o.FeasibilityTol, ...
        'IntFeasTol', o.IntFeasTol, 'TimeLimit', o.TimeLimit, ...
        'OutputFlag', double(o.Verbose), 'NumericFocus', o.NumericFocus, ...
        'ScaleFlag', o.ScaleFlag);
    if o.LegacyMethod, params.Method = 5; end

    wt = solveLinear(Aeq, b, model.lb, model.ub, unitObjective(n,bio,-1), params, 'maximum biomass');
    optBio = wt.x(bio);
    product = solveLinear(Aeq, b, model.lb, model.ub, unitObjective(n,target,-1), params, 'maximum product');
    refLb = model.lb; refLb(bio) = max(refLb(bio), o.ReferenceGrowthFraction * optBio);
    productAtGrowth = solveLinear(Aeq, b, refLb, model.ub, unitObjective(n,target,-1), params, 'product at reference growth');
    maxProductAtGrowth = productAtGrowth.x(target);
    if maxProductAtGrowth <= o.ZeroFluxTolerance
        error('okoSolve:NoReferenceProduction', ...
            'Target flux at reference growth is zero; OKO cannot construct a production reference.');
    end

    enzymeLb = model.lb;
    enzymeLb(bio) = max(enzymeLb(bio), o.EnzymeReferenceGrowthFraction * optBio);
    enzymeLb(target) = max(enzymeLb(target), o.ReferenceProductFraction * maxProductAtGrowth);
    enzymeObj = zeros(n,1); enzymeObj(prepared.enzymeVarIdx) = prepared.enzymeVarSign;
    enzymeRef = solveLinear(Aeq, b, enzymeLb, model.ub, enzymeObj, params, 'minimum enzyme usage');
    eLevel = prepared.enzymeVarSign .* enzymeRef.x(prepared.enzymeVarIdx);
    minE = max((1-o.EnzymeTolerance) .* eLevel, 0);
    maxE = (1+o.EnzymeTolerance) .* eLevel;

    pairMask = true(height(prepared.pairs),1);
    if isPlus, pairMask = prepared.pairs.hasInterval; end
    pairs = prepared.pairs(pairMask,:);
    if isempty(pairs)
        error('okoSolve:NoBoundedPairs', 'No model kcat pairs matched the supplied OKO+ intervals.');
    end
    groups = remapGroups(prepared.complexGroups, find(pairMask));
    problem = buildMilp(Aeq, b, model.lb, model.ub, prepared.enzymeVarIdx, prepared.enzymeVarSign, ...
        minE, maxE, pairs, groups, isPlus, o);

    warmLb = problem.lb; warmUb = problem.ub;
    warmLb(bio) = max(warmLb(bio), o.EnzymeReferenceGrowthFraction * optBio);
    warmUb(target) = min(warmUb(target), o.WarmupProductFold * maxProductAtGrowth);
    warmUb(bio) = min(warmUb(bio), o.WarmupGrowthUpperFraction * optBio);
    warmObj = zeros(numel(problem.lb),1); warmObj(target) = -1;
    if contains(prepared.profile,'legacy') && ~isPlus
        fake=zeros(numel(problem.lb),1); baseStart=enzymeRef.x;
        baseStart(productAtGrowth.x(1:n)<0)=0; fake(1:n)=baseStart;
        fake(problem.index.yK)=1; problem.start=fake;
    end
    warm = solveGurobi(problem, warmLb, warmUb, warmObj, params, 'OKO warm-up');
    if strcmpi(warm.status, 'INFEASIBLE') || strcmpi(warm.status, 'INF_OR_UNBD')
        warmLb(bio) = max(model.lb(bio), o.RelaxedGrowthFraction * optBio);
        warmUb(bio) = problem.ub(bio);
        warm = solveGurobi(problem, warmLb, warmUb, warmObj, params, 'relaxed OKO warm-up');
    end
    assertSolution(warm, 'OKO warm-up');

    finalObj = zeros(numel(problem.lb),1);
    finalObj(problem.index.yK) = -1; % maximize unchanged kcats
    if isPlus, finalObj(problem.index.yE) = -o.AbundanceWeight; end % and unchanged abundance (w=10 per paper)
    final = struct('status','INFEASIBLE');
    productFraction = o.EngineeredProductFraction;
    growthFraction = o.EngineeredGrowthFraction;
    while productFraction >= o.MinimumEngineeredFraction-1e-12
        finalLb = warmLb; finalUb = warmUb;
        finalLb(target) = max(finalLb(target), productFraction * warm.x(target));
        finalLb(bio) = max(finalLb(bio), growthFraction * warm.x(bio));
        finalStart=warm.x;
        if contains(prepared.profile,'legacy')
            negativeBase=find(finalStart(1:n)<0); finalStart(negativeBase)=0;
        end
        problem.start = clip(finalStart, finalLb, finalUb);
        final = solveGurobi(problem, finalLb, finalUb, finalObj, params, 'OKO intervention minimization');
        if isfield(final,'x') && ~any(strcmpi(final.status,{'INFEASIBLE','INF_OR_UNBD'})), break; end
        productFraction = productFraction-o.EngineeredFractionStep;
        growthFraction = growthFraction-o.EngineeredFractionStep;
    end
    assertSolution(final, 'OKO intervention minimization');
    o.ActualEngineeredProductFraction = productFraction;
    o.ActualEngineeredGrowthFraction = growthFraction;

    result = formatResult(prepared, pairs, pairMask, mode, o, wt, product, ...
        productAtGrowth, enzymeRef, warm, final, problem, minE, maxE);
    if ~isempty(o.OutputFile)
        outputFile = char(o.OutputFile);
        parent = fileparts(outputFile); if ~isempty(parent) && ~exist(parent,'dir'), mkdir(parent); end
        writetable(result.changes, outputFile);
        result.outputFile = outputFile;
    end
    if ~isempty(o.MatFile)
        matFile = char(o.MatFile);
        parent = fileparts(matFile); if ~isempty(parent) && ~exist(parent,'dir'), mkdir(parent); end
        result.matFile = matFile; save(matFile, 'result');
    end
end

function problem = buildMilp(Aeq, b, lb, ub, enzymeVars, enzymeSigns, minE, maxE, pairs, groups, isPlus, o)
    n = size(Aeq,2); p = height(pairs); ne = numel(enzymeVars);
    yEcount = ne * isPlus;
    idx.yE = n + (1:yEcount);
    idx.dp = n + yEcount + (1:p);
    idx.yK = n + yEcount + p + (1:p);
    total = n + yEcount + 2*p;

    E = sparse(pairs.metIdx, 1:p, 1, size(Aeq,1), p);
    AeqFinal = [Aeq, sparse(size(Aeq,1), yEcount), E, sparse(size(Aeq,1),p)];
    beqFinal = b;
    for g = 1:numel(groups)
        members = groups{g};
        for j = 1:numel(members)-1
            row = sparse(1,total); row(idx.dp(members(j)))=1; row(idx.dp(members(j+1)))=-1;
            AeqFinal(end+1,:) = row; beqFinal(end+1,1) = 0; %#ok<AGROW>
        end
    end

    rows = {}; rhs = {};
    if isPlus
        R1 = sparse(ne,total); R2 = sparse(ne,total);
        R1(sub2ind(size(R1),(1:ne)',enzymeVars(:))) = -enzymeSigns;
        R2(sub2ind(size(R2),(1:ne)',enzymeVars(:))) = enzymeSigns;
        R1(sub2ind(size(R1),(1:ne)',idx.yE(:))) = o.BigM;
        R2(sub2ind(size(R2),(1:ne)',idx.yE(:))) = o.BigM;
        rows{end+1} = R1; rhs{end+1} = o.BigM*ones(ne,1)-minE;
        rows{end+1} = R2; rhs{end+1} = o.BigM*ones(ne,1)+maxE;
    else
        positive = enzymeSigns > 0; negative = ~positive;
        lb(enzymeVars(positive)) = max(lb(enzymeVars(positive)), minE(positive));
        ub(enzymeVars(positive)) = min(ub(enzymeVars(positive)), maxE(positive));
        lb(enzymeVars(negative)) = max(lb(enzymeVars(negative)), -maxE(negative));
        ub(enzymeVars(negative)) = min(ub(enzymeVars(negative)), -minE(negative));
    end

    coeff = pairs.coefficient;
    R1 = sparse(p,total); R2 = sparse(p,total);
    R1(sub2ind(size(R1),(1:p)',pairs.rxnIdx)) = -o.Significance/(1-o.Significance);
    R2(sub2ind(size(R2),(1:p)',pairs.rxnIdx)) = -o.Significance/(1+o.Significance);
    R1(sub2ind(size(R1),(1:p)',idx.dp(:))) = 1./coeff;
    R2(sub2ind(size(R2),(1:p)',idx.dp(:))) = -1./coeff;
    R1(sub2ind(size(R1),(1:p)',idx.yK(:))) = o.BigM;
    R2(sub2ind(size(R2),(1:p)',idx.yK(:))) = o.BigM;
    rows(end+1:end+2) = {R1,R2}; rhs(end+1:end+2) = {o.BigM*ones(p,1),o.BigM*ones(p,1)};

    R3 = sparse(p,total);
    R3(sub2ind(size(R3),(1:p)',pairs.rxnIdx)) = o.PositiveGuard-1;
    R3(sub2ind(size(R3),(1:p)',idx.dp(:))) = -1./coeff;
    rows{end+1}=R3; rhs{end+1}=zeros(p,1);
    R4=sparse(p,total); R5=sparse(p,total);
    R4(sub2ind(size(R4),(1:p)',pairs.rxnIdx))=-o.BigM;
    R5(sub2ind(size(R5),(1:p)',pairs.rxnIdx))=-o.BigM;
    R4(sub2ind(size(R4),(1:p)',idx.dp(:)))=-1;
    R5(sub2ind(size(R5),(1:p)',idx.dp(:)))=1;
    rows(end+1:end+2)={R4,R5}; rhs(end+1:end+2)={o.ZeroFluxGuard*ones(p,1),o.ZeroFluxGuard*ones(p,1)};

    R6=sparse(p,total); R7=sparse(p,total);
    if isPlus
        zMin = pairs.kcatMin ./ pairs.factor; zMax = pairs.kcatMax ./ pairs.factor;
        R6(sub2ind(size(R6),(1:p)',pairs.rxnIdx)) = zMin + 1./coeff;
        R7(sub2ind(size(R7),(1:p)',pairs.rxnIdx)) = -zMax - 1./coeff;
        R6(sub2ind(size(R6),(1:p)',idx.dp(:))) = (1./coeff).*zMin;
        R7(sub2ind(size(R7),(1:p)',idx.dp(:))) = (-1./coeff).*zMax;
    else
        R6(sub2ind(size(R6),(1:p)',pairs.rxnIdx)) = 1-o.FoldLimit;
        R7(sub2ind(size(R7),(1:p)',pairs.rxnIdx)) = 1/o.FoldLimit-1;
        R6(sub2ind(size(R6),(1:p)',idx.dp(:))) = 1./coeff;
        R7(sub2ind(size(R7),(1:p)',idx.dp(:))) = -1./coeff;
    end
    rows(end+1:end+2)={R6,R7}; rhs(end+1:end+2)={zeros(p,1),zeros(p,1)};
    Aineq = vertcat(rows{:}); bineq = vertcat(rhs{:});

    problem = struct('A',[Aineq;AeqFinal], 'rhs',[bineq;beqFinal], ...
        'sense',[repmat('<',size(Aineq,1),1);repmat('=',size(AeqFinal,1),1)], ...
        'lb',[lb;zeros(yEcount,1);-o.DeltaBound*ones(p,1);zeros(p,1)], ...
        'ub',[ub;ones(yEcount,1);o.DeltaBound*ones(p,1);ones(p,1)], ...
        'vtype',repmat('C',total,1),'index',idx);
    problem.vtype(idx.yK)='B'; if isPlus, problem.vtype(idx.yE)='B'; end
end

function result = formatResult(prepared, pairs, pairMask, mode, o, wt, product, prodRef, enzRef, warm, final, problem, minE, maxE)
    flux = final.x(pairs.rxnIdx); delta = final.x(problem.index.dp);
    deltaCoeff = nan(height(pairs),1); active = abs(flux)>o.ZeroFluxTolerance;
    deltaCoeff(active) = delta(active)./flux(active);
    newKcat = pairs.factor ./ (-(pairs.coefficient + deltaCoeff));
    yK = round(final.x(problem.index.yK));
    changed = yK==0 & active & isfinite(newKcat) & newKcat>0 & ...
        abs(log(newKcat./pairs.kcat)) > o.ReportLogTolerance;
    changes = pairs(changed,{'pairID','enzymeID','rxnID','rxnName','kcat','kcatMin','kcatMax'});
    changes.newKcat = newKcat(changed);
    changes.foldChange = newKcat(changed)./pairs.kcat(changed);
    changes.logFold = log(changes.foldChange);
    abundanceChanges = table();
    if strcmp(mode,'oko-plus')
        yE=round(final.x(problem.index.yE));
        e=prepared.enzymeVarSign.*final.x(prepared.enzymeVarIdx);
        ch=yE==0;
        abundanceChanges=table(prepared.enzymeVarIdx(ch),e(ch),minE(ch),maxE(ch), ...
            'VariableNames',{'reactionIndex','engineeredFlux','referenceMin','referenceMax'});
    end
    result=struct('algorithm',mode,'profile',prepared.profile,'config',o, ...
        'biomassRxn',prepared.biomassRxn,'targetRxn',prepared.targetRxn, ...
        'diagnostics',struct('wtBiomass',wt.x(prepared.biomassIdx), ...
        'maxProduct',product.x(prepared.targetIdx),'referenceProduct',prodRef.x(prepared.targetIdx), ...
        'referenceEnzymeObjective',enzRef.objval,'engineeredBiomass',final.x(prepared.biomassIdx), ...
        'engineeredProduct',final.x(prepared.targetIdx),'warmupProduct',warm.x(prepared.targetIdx), ...
        'solverStatus',final.status,'eligiblePairs',height(pairs),'totalPairs',height(prepared.pairs), ...
        'intervalMatchedPairs',nnz(pairMask)), ...
        'changes',changes,'abundanceChanges',abundanceChanges,'pairs',pairs, ...
        'outputFile','','matFile','');
end

function sol=solveLinear(A,b,lb,ub,obj,params,label)
    m=struct('A',A,'rhs',full(b),'sense',repmat('=',size(A,1),1),'lb',lb,'ub',ub, ...
        'obj',obj,'vtype',repmat('C',size(A,2),1));
    sol=gurobi(m,params); assertSolution(sol,label);
end
function sol=solveGurobi(problem,lb,ub,obj,params,label)
    m=problem; m.lb=lb; m.ub=ub; m.obj=obj; m=rmfield(m,'index'); sol=gurobi(m,params);
    if ~isfield(sol,'status'), error('okoSolve:SolverFailure','Gurobi returned no status for %s.',label); end
end
function assertSolution(sol,label)
    good={'OPTIMAL','SUBOPTIMAL','TIME_LIMIT'};
    if ~isfield(sol,'status') || ~any(strcmpi(sol.status,good)) || ~isfield(sol,'x')
        status='UNKNOWN'; if isfield(sol,'status'),status=sol.status;end
        error('okoSolve:SolverFailure','%s failed with status %s.',label,status);
    end
end
function f=unitObjective(n,idx,value), f=zeros(n,1); f(idx)=value; end
function x=clip(x,lb,ub), x=max(x,lb); x=min(x,ub); end
function groups=remapGroups(oldGroups,selected)
    groups={}; map=zeros(max([selected(:);0]),1); map(selected)=1:numel(selected);
    for i=1:numel(oldGroups)
        old=oldGroups{i}; old=old(old<=numel(map)); members=map(old); members=members(members>0);
        if numel(members)>1, groups{end+1,1}=members(:)'; end %#ok<AGROW>
    end
end
function o=defaults(in,profile)
    legacy=contains(lower(profile),'legacy');
    d=struct('FluxScale',1,'NormalizeRows',true,'LegacyMethod',legacy, ...
        'ReferenceGrowthFraction',0.99,'EnzymeReferenceGrowthFraction',0.98, ...
        'ReferenceProductFraction',0.99,'EnzymeTolerance',0.1,'Significance',1e-8, ...
        'FoldLimit',10,'BigM',1e6,'DeltaBound',1e12,'PositiveGuard',1e-12, ...
        'ZeroFluxGuard',1e-6,'ZeroFluxTolerance',1e-9,'WarmupProductFold',2, ...
        'WarmupGrowthUpperFraction',0.99,'RelaxedGrowthFraction',0.5, ...
        'EngineeredProductFraction',0.99,'EngineeredGrowthFraction',0.99, ...
        'MinimumEngineeredFraction',0.90,'EngineeredFractionStep',0.01, ...
        'AbundanceWeight',10, ...   % w in max(sum y_k + w*sum y_e); paper Methods p.15 sets w=10
        'ReportLogTolerance',1e-6,'FeasibilityTol',1e-6,'IntFeasTol',1e-6, ...
        'NumericFocus',0,'ScaleFlag',-1,'TimeLimit',900,'Verbose',false, ...
        'OutputFile','','MatFile','');
    if legacy, d.FluxScale=1000; end
    o=d; names=fieldnames(in); for i=1:numel(names), o.(names{i})=in.(names{i}); end
end
