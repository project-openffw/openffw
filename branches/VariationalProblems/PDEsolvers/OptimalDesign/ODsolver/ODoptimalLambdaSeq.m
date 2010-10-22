function p = ODoptimalLambdaSeq(p)
% author: David Guenther 
% compute the problem specific optimal lambda

%% create fine mesh
nrRefine = 1;

p.level(1).level = 1;
postProc = p.statics.postProc;
enumerate = p.statics.enumerate;
mark = p.statics.mark;
refine = p.statics.refine;
pdeSolver = p.params.pdeSolver;
prolong = str2func([pdeSolver,'prolong']);
getNonLinearSolution = str2func([pdeSolver,'getNonLinearSolution']);

p.problem.lambda = 0.1;
p.problem.t1 = 1;
p.problem.t2 = 2;

for k = 1:nrRefine
    p = enumerate(p);
    p = mark(p);
    p = refine(p);
    p = enumerate(p);
    p = prolong(p,k+1);
    p = postProc(p);
    p.level(end).level = k+1;
end
%0.01315561749642   0.0168
lambdaVals = [0.013,0.0131,0.0132,0.014,0.015,0.016,0.0168,0.017];
%lambdaVals = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
fprintf('\n');

mu1 = p.problem.mu1;
mu2 = p.problem.mu2;

format long

options = p.params.options;

degree = p.params.rhsIntegtrateExactDegree;
output = [0 0];
%% compute optimal lambda
for lambda = lambdaVals
    t1 = (2*lambda*mu1/mu2).^(1/2);
    t2 = mu2/mu1*t1;
    p.problem.lambda = lambda;
    p.problem.t1 = t1;
    p.problem.t2 = t2;
    
    n4e = p.level(end).geom.n4e;
%     p.level(nrRefine+1).f4e = integrate(n4e,nrRefine+1,degree,@funcHandleRHSVolume,p);
%    p.level(nrRefine+1).f4e = integrate(n4e,nrRefine+1,degree,@RHS,p);
    p.level(nrRefine+1).f4e = integrate(n4e,nrRefine+1,degree,@RHS2,p);
    freeNodes = p.level(end).enum.freeNodes;
    x0 = p.level(end).x;
    p.level(end).x0 = x0;
    % find x s.t. E(x) = 0 with E given in getFuncVal.m
    [x,fval,exitflag,outputSolve,jacobian] = fsolve(getNonLinearSolution,x0(freeNodes),options,p);
    x0(freeNodes) = x;

    p.level(end).x = x0;
    p.level(end).fval = fval;
    p.level(end).jacobian = jacobian;
    p.level(end).exitflag = exitflag;
    p.level(end).output = outputSolve;

    p = postProc(p);

    k = 1;
    energy = getEnergy(p);

    output = [output;[lambda energy]];
end

output

format long
fprintf('\nOptimal Lambda for this problem is:\n lambda = % 3.4f',lambda);
fprintf('\nwith energy:\n energy = % 3.4f\n',energy);
format short
%% supply the discrete energy
function val = getEnergy(p)

energy_h = p.statics.energy_h;
lvl = size(p.level,2);
n4e = p.level(lvl).geom.n4e;

intVal = integrate(n4e,lvl,10,energy_h,p);
val = sum(intVal);

function val = RHS(x,y,curElem,lvl,p)

sigma0 = p.problem.sigma0;
stressBasis = p.level(lvl).enum.grad4e;

evalSigma = sigma0(x,y,curElem,lvl,p);
evalBasis = stressBasis(:,:,curElem);

val = evalBasis*evalSigma';
val = reshape(val,[3 1 length(x)]);

function val = RHS2(x,y,curElem,lvl,p)

f = p.problem.f;
basis = p.statics.basisU;

evalBasis = basis(x,y,curElem,lvl,p);
nrBasis = size(evalBasis,2);

evalF = f(x,y,curElem,lvl,p)*ones(1,nrBasis);

val = reshape(evalF'.*evalBasis',[nrBasis 1 length(x)]);
