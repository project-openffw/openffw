function p = TWdirectFsolve(p)

% author: David Guenther, Lena Noack 
% load abort criteria
maxLevel = loadField('p.params','maxLevel',p,100);
maxNrDoF = loadField('p.params','maxNrDoF',p,500);
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;

DWRbool = loadField('p.params','DWRbool',p,0);


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

if DWRbool == 1
      GOAL = p.params.GOAL;
      %GOAL_IM: integral middle:   J(tau,v) = 1/|Omega1|  int_Omega1 v dx
      %GOAL_DP: dual problem:      J(tau,v) = 1/|Omega1|  int_Omega1 Dv dx
      %GOAL_PR: for rel-eff proof: J(tau,v) = int_Omega tau(sigma-sigma_ell) dx
      if strcmp(GOAL,'GOAL_IM')
        getNonLinearSolutionDWR = str2func([pdeSolver,'getNonLinearSolutionDWRIM']);
      elseif strcmp(GOAL,'GOAL_DP')
        getNonLinearSolutionDWR = str2func([pdeSolver,'getNonLinearSolutionDWRDP']);
      else
        getNonLinearSolutionDWR = str2func([pdeSolver,'getNonLinearSolutionDWRPR']);
      end
end
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

    [x,fval,exitflag,output,jacobian] = fsolve(getNonLinearSolution,x0(freeNodes),options,p);
    x0(freeNodes) = x;

    p.level(end).x = x0;
    p.level(end).fval = fval;
    p.level(end).jacobian = jacobian;
    p.level(end).exitflag = exitflag;
    p.level(end).output = output;

    if DWRbool == 1
        DWRx0 = p.level(end).DWRx;
        p.level(end).DWRx0 = DWRx0;
        DWRf4e = integrate(n4e,curLevel+1,degree,@fDual,p);
        p.level(curLevel+1).DWRf4e = DWRf4e;
        [DWRx] = fsolve(getNonLinearSolutionDWR,DWRx0(freeNodes),options,p);
        DWRx0(freeNodes) = DWRx;
%        fprintf('\n level = %d, |x| = %.4g \t, |DWRx| = %.4g',curLevel,norm(x0,2),norm(DWRx0,2));

        p.level(end).DWRx = DWRx0;
    end

    p = postProc(p);
    p = estimate(p);
    curLevel = curLevel + 1;
    p.level(end).level = curLevel;
    nrDoF = p.level(end).nrDoF;
    estimatedError = p.level(end).estimatedError;
    fprintf('\n Level = %d \t Error = %.2g \t DoF = %d\n',curLevel,estimatedError,nrDoF);
end

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

function val = fDual(x,y,curElem,lvl,p)

basis = p.statics.basisU;

evalBasis = basis(x,y,curElem,lvl,p);
nrBasis = size(evalBasis,2);

evalF = zeros(size(x,1),nrBasis);

GOAL = p.params.GOAL;
if strcmp(GOAL,'GOAL_IM') %mean integral of v

%f=1 on Omega2 (Square domain) and 0 on O1:
%  _____
% |  |  |
% |O1|O2|
% |__|__|
%
index1 = find( x >= 0.5 );
evalF(index1,:) = 4/3;

elseif strcmp(GOAL,'GOAL_DP') %mean integral of Dv

%f=1 on Omega2 (Square domain) and 0 on O1:
%  ___________
% |  |        |
% |  |        |
% |O1|   O2   |
% |  |        |
% |__|________|
%

index1 = find( x <= 0.25 );
evalF(index1,:) = 8/3;


else %GOAL_PR

%f=1 on Omega2 (Square domain) and 0 on O1:
%  _____
% |     |
% |  O  |
% |__ __|
%
index1 = find( x >= 0 );
evalF(index1,:) = 1;%2/3;% 1 inseat 2/3, since no mean integral was calculated

end

val = reshape(evalF'.*evalBasis',[nrBasis 1 length(x)]);