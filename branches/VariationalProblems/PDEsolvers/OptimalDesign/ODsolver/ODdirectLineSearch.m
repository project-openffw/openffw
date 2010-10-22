function p = ODdirectLineSearch(p)
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
while(nrDoF < maxNrDoF && curLevel <= maxLevel)
    p = mark(p);
    p = refine(p);
    p = enumerate(p);
    p = prolong(p,curLevel+1);

    n4e = p.level(end).geom.n4e;
    p.level(curLevel+1).f4e = integrate(n4e,curLevel+1,degree,@funcHandleRHSVolume,p);
    freeNodes = p.level(end).enum.freeNodes;
    x0 = p.level(end).x;
    p.level(end).x0 = x0;
    % find x s.t. E(x) = 0 with E given in getFuncVal.m
    [x,fval,jacobian,t] = ODlineSearch(x0(freeNodes),curLevel,p);
    x0(freeNodes) = x;

    p.level(end).x = x0;
    p.level(end).fval = fval;
    p.level(end).jacobian = jacobian;

    p = postProc(p);
    p = estimate(p);
    curLevel = curLevel + 1;
    p.level(end).level = curLevel;
    nrDoF = p.level(end).nrDoF;
    estimatedError = p.level(end).estimatedError;
    fprintf('\n Level = %d \t Error = %.2g \t DoF = %d\n',curLevel+1,estimatedError,nrDoF);
end
