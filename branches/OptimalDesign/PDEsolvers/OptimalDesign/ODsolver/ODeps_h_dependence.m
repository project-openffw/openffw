function p = ODeps_h_dependence(p)
% author: David Guenther 
% Copyright 2007 David Guenther
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%

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

    if(curLevel == 1)
        p.level(1).markedEdges = true( length(p.level(1).markedEdges),1 );
    end
    
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
function val = RHS1(points,curElem,lvl,p)

sigma0 = p.problem.sigma0;
stressBasis = p.statics.stressBasis;

evalSigma = sigma0(points,curElem,lvl,p);
evalBasis = stressBasis(points,curElem,lvl,p);

val = matMul(evalBasis,reshape(evalSigma',[2 1 length(points(:,1))]));

%% compute RHS:= f w
function val = RHS2(points,curElem,lvl,p)

f = p.problem.f;
basis = p.statics.basisU;

evalBasis = basis(points,curElem,lvl,p);
nrBasis = size(evalBasis,2);

evalF = f(points,curElem,lvl,p)*ones(1,nrBasis);

val = reshape(evalF'.*evalBasis',[nrBasis 1 length(points(:,1))]);
