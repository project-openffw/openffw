function [val,jacobian] = ODP1getNonLinearSolution(x,p)
%author: David Guenther

%% integrate Dirichlet boundary in x and save the variable x
x0 = p.level(end).x0;
freeNodes = p.level(end).enum.freeNodes;
x0(freeNodes) = x;
p.level(end).x0 = x0;
p.level(end).x = x0;

postProc = p.statics.postProc;
p = postProc(p);

pdeSolver = p.params.pdeSolver;
%% get function-value E(x)
getFuncVal = str2func([pdeSolver,'getFuncVal']);
p = getFuncVal(p);
val = p.level(end).funcVal;

%% get the jacobian of E(x)
if nargout > 1
    getJacobian = str2func([pdeSolver,'getJacobian']);
    p = getJacobian(p);
    jacobian = p.level(end).jacobi;
end
