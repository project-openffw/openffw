function [fval,gradient,hessian] = ODRTgetNL4fminunc(x,p)
% author: David Guenther 

%% integrate Dirichlet boundary in x and save the variable x
p.level(end).x0 = x;
p.level(end).x = x;

postProc = p.statics.postProc;
p = postProc(p);

pdeSolver = p.params.pdeSolver;
%% get discrete energy
n4e = p.level(end).geom.n4e;
% energy_h = p.statics.energy_h;
lvl = length(p.level);
energy4e = integrate(n4e,lvl,10,@getConjFunctional,p);
% energy4e = integrate(n4e,lvl,10,energy_h,p);
fval = sum(energy4e);

%% get function-value E(x)
if nargout > 1
    getFuncVal = str2func([pdeSolver,'getFuncVal']);
    p = getFuncVal(p);
    gradient = p.level(end).funcVal;
end
%% get the jacobian of E(x)
if nargout > 2
    getJacobian = str2func([pdeSolver,'getJacobian']);
    p = getJacobian(p);
    hessian = p.level(end).jacobi;
end


function val = getConjFunctional(points,curElem,lvl,p)

nonLinearConjExact = p.problem.nonLinearConjExact;
sigma_h = p.statics.sigma_h;

evalSigma = sigma_h(points,curElem,lvl,p);
absSigma = (evalSigma(:,1).^2 + evalSigma(:,2).^2).^(1/2);
evalNonLinear = nonLinearConjExact(absSigma,curElem,lvl,p);

val = reshape(evalNonLinear,[1 1 length(points(:,1))]);