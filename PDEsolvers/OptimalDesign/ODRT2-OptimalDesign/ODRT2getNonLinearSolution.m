function [val,jacobian] = ODRTgetNonLinearSolution(x,p)
% author: David Guenther 

%% save the variable x
p.level(end).x0 = x;
p.level(end).x = x;

postProc = p.statics.postProc;
p = postProc(p);

pdeSolver = p.params.pdeSolver;
%% get function-value E(x)
getFuncVal = str2func(['ODRTgetFuncVal']);
%getFuncVal = str2func([pdeSolver,'getFuncVal']);
p = getFuncVal(p);
val = p.level(end).funcVal;

%% get the jacobian of E(x)
if nargout > 1
    getJacobian = str2func(['ODRTgetJacobian']);
%    getJacobian = str2func([pdeSolver,'getJacobian']);
    p = getJacobian(p);
    jacobian = p.level(end).jacobi;
end
