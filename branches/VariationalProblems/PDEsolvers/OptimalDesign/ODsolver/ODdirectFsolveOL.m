function p = ODdirectFsolve(p)
% author: David Guenther 
% load abort criteria
maxLevel = loadField('p.params','maxLevel',p,100);
maxNrDoF = loadField('p.params','maxNrDoF',p,500);

nrDoF = p.level(1).nrDoF;
curLevel = 1;
p.level(1).level = curLevel;

enumerate = p.statics.enumerate;
mark = p.statics.mark;
refine = p.statics.refine;
postProc = p.statics.postProc;
estimate = p.statics.estimate;

options = p.params.options;
pdeSolver = p.params.pdeSolver;
prolong = str2func([pdeSolver,'prolong']);
getNonLinearSolution = str2func([pdeSolver,'getNonLinearSolution']);

degree = p.params.rhsIntegtrateExactDegree;
RHS = p.params.RHS;
while(nrDoF < maxNrDoF && curLevel <= maxLevel)
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

    warning off
    [x,fval,exitflag,output,jacobian] = fsolve(getNonLinearSolution,x0(freeNodes),options,p);
    warning on
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
%    fprintf('\n Level = %d \t Error = %.2g \t DoF = %d\n',curLevel+1,estimatedError,nrDoF);

end
    energy1 = getEnergy(p);
    lambda = p.problem.lambda;
    mu1 = p.problem.mu1;
    mu2 = p.problem.mu2;
    dP = loadField('p.params','distributeParam',p,0.5);
    fprintf('Lambda = %.5g \t Energy = %.5g \t DoF = %d \t mu1 = %d \t mu2 = %d \t xi = %.2g\n',...
        lambda,energy1,nrDoF,mu1,mu2,dP);
%    fprintf('\n Level = %d \t Lambda = %.5g \t Energy = %.5g \t DoF = %d\n',curLevel+1,lambda,energy1,nrDoF);

function val = RHS1(x,y,curElem,lvl,p)

sigma0 = p.problem.sigma0;
stressBasis = p.statics.stressBasis;

evalSigma = sigma0(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

val = matMul(evalBasis,reshape(evalSigma',[2 1 length(x)]));

function val = RHS2(x,y,curElem,lvl,p)

f = p.problem.f;
basis = p.statics.basisU;

evalBasis = basis(x,y,curElem,lvl,p);
nrBasis = size(evalBasis,2);

evalF = f(x,y,curElem,lvl,p)*ones(1,nrBasis);

val = reshape(evalF'.*evalBasis',[nrBasis 1 length(x)]);



%% supply the discrete energy
function val = getEnergy(p)

energy_h = p.statics.energy_h;
lvl = size(p.level,2);
n4e = p.level(lvl).geom.n4e;

intVal = integrate(n4e,lvl,10,energy_h,p);
val = sum(intVal);
