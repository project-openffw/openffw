function [fval,gradient,hessian] = P1getNL4fminunc(x,p)
%author: David Guenther, Lena Noack

%% integrate Dirichlet boundary in x and save the variable x
x0 = p.level(end).x0;
freeNodes = p.level(end).enum.freeNodes;
x0(freeNodes) = x;
p.level(end).x0 = x0;
p.level(end).x = x0;

postProc = p.statics.postProc;
p = postProc(p);

pdeSolver = p.params.pdeSolver;
%% get discrete energy
n4e = p.level(end).geom.n4e;
energy_h = p.statics.energy_h;
lvl = length(p.level);
energy4e = integrate(n4e,lvl,10,energy_h,p);
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
