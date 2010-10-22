function [val,jacobian] = ODRTDWRgetNonLinearSolution(x,p)
% author: David Guenther, Lena Noack

%% save the variable x
p.level(end).DWRx0 = x;
p.level(end).DWRx = x;

postProc = p.statics.postProc;
p = postProc(p);

pdeSolver = p.params.pdeSolver;
%% get function-value E(x)
getFuncVal = str2func(['ODRTDWRgetFuncValIM']);
p = getFuncVal(p);
val = p.level(end).funcVal;

%% get the jacobian of E(x)
if nargout > 1
    getJacobian = str2func(['ODRTDWRgetJacobianIM']);
    p = getJacobian(p);
    jacobian = p.level(end).jacobi;
end
