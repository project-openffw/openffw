function p = ODeps_h_dependence(p)
% author: David Guenther 

% epsilon is mesh dependent

maxLevel = loadField('p.params','maxLevel',p,100);
maxNrDoF = loadField('p.params','maxNrDoF',p,500);

nrDoF = p.level(1).nrDoF;
curLevel = 1;
p.level(1).level = curLevel;
k = 1;

fprintf('\n');

options = p.params.options;

enumerate = p.statics.enumerate;
mark = p.statics.mark;
refine = p.statics.refine;
postProc = p.statics.postProc;
estimate = p.statics.estimate;

pdeSolver = p.params.pdeSolver;
prolong = str2func([pdeSolver,'prolong']);
getNonLinearSolution = str2func([pdeSolver,'getNonLinearSolution']);

degree = p.params.rhsIntegtrateExactDegree;
RHS = p.params.RHS;

fprintf('\n');

while( (nrDoF < maxNrDoF && curLevel <= maxLevel) )
    p = mark(p);
    p = refine(p);
    p = enumerate(p);
    p = prolong(p,curLevel+1);
    n4e = p.level(end).geom.n4e;
    if strcmp(RHS,'RHS1')
        p.level(curLevel+1).f4e = integrate(n4e,curLevel+1,degree,@RHS1,p);
    else
        p.level(curLevel+1).f4e = integrate(n4e,curLevel+1,degree,@RHS2,p);
    end

    freeNodes = p.level(end).enum.freeNodes;

    x0 = p.level(end).x;
    p.level(end).x0 = x0;
    % find x s.t. E(x) = 0 with E given in getFuncVal.m
    [x,fval,exitflag,output,jacobian] = fsolve(getNonLinearSolution,x0(freeNodes),options,p);
    x0(freeNodes) = x;

    p.level(end).x = x0;
    p.level(end).fval = fval;
    p.level(end).jacobian = jacobian;
    p.level(end).exitflag = exitflag;
    p.level(end).output = output;

    p = postProc(p);

    p = estimate(p);
    curLevel = curLevel + 1;
    p.level(end).level = curLevel;
    nrDoF = p.level(end).nrDoF;
    estimatedError = p.level(end).estimatedError;
    fprintf('\n Level = %d \t Error = %.2g \t DoF = %d\n',curLevel+1,estimatedError,nrDoF);
end

%% compute RHS:= f \nabla w
function val = RHS1(x,y,curElem,lvl,p)

sigma0 = p.problem.sigma0;
stressBasis = p.statics.stressBasis;

evalSigma = sigma0(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

val = matMul(evalBasis,reshape(evalSigma',[2 1 length(x)]));

%% compute RHS:= f w
function val = RHS2(x,y,curElem,lvl,p)

f = p.problem.f;
basis = p.statics.basisU;

evalBasis = basis(x,y,curElem,lvl,p);
nrBasis = size(evalBasis,2);

evalF = f(x,y,curElem,lvl,p)*ones(1,nrBasis);

val = reshape(evalF'.*evalBasis',[nrBasis 1 length(x)]);
