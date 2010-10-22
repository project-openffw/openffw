function p = directFsolve(p)
% author: David Guenther 
% load abort criteria
maxLevel = loadField('p.params','maxLevel',p,100);
maxNrDoF = loadField('p.params','maxNrDoF',p,500);
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
    fprintf('\n Level = %d \t Error = %.2g \t DoF = %d\n',curLevel+1,estimatedError,nrDoF);
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
problem = p.params.problem.name; 
if strcmp(GOAL,'GOAL_IM') %mean integral of v

if (strcmp(problem,'AbsValue_Square') || strcmp(problem,'AbsValue_Square_exact'))
%f=1 on Omega2 (Square domain) and 0 on O1:
%  _____
% |  |  |
% |O1|O2|
% |__|__|
%
index1 = find( x >= 0 );
evalF(index1,:) = 0.5;
else
%f=1 on Omega2 (Lshape domain) and 0 on O1,O3:
%  _________
% |    |    |
% | O1 | O2 |
% |____|____|
% |    |
% | O3 |
% |____|
%
index1 = find( x >= 0 );
evalF(index1,:) = 1;
end

elseif strcmp(GOAL,'GOAL_DP') %mean integral of Dv

%f=1 on Omega2 (Square domain) and 0 on O1:
%  ___________
% |  |        |
% |  |        |
% |O1|   O2   |
% |  |        |
% |__|________|
%

% f=1 on Omega 1 (Lshape domain) and 0 else
%  ___________
% |  |        |
% |  |   O2   |
% |O1|   _____|
% |  |  |
% |  |  |
% |__|__|
%

index1 = find( x <= -0.5 );
evalF(index1,:) = 1;


else %GOAL_PR

if (strcmp(problem,'AbsValue_Square') || strcmp(problem,'AbsValue_Square_exact'))
%f=1 on Omega2 (Square domain) and 0 on O1:
%  _____
% |     |
% |  O  |
% |__ __|
%
index1 = find( x >= 0 );
evalF(index1,:) = 1;%0.25;
else
%f=1 on Omega2 (Lshape domain) and 0 on O1,O3:
%  _________
% |         |
% |  O      |
% |     ____|
% |    |
% |    |
% |____|
%
index1 = find( x >= 0 );
evalF(index1,:) = 1;%1/3;
end

end

val = reshape(evalF'.*evalBasis',[nrBasis 1 length(x)]);
