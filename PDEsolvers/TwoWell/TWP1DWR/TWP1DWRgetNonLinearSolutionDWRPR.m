function [val,jacobian] = P1getNonLinearSolutionDWRPR(DWRx,p)
%author: David Guenther, Lena Noack

%% integrate Dirichlet boundary in x and save the variable x
%x0 = p.level(end).x0;
%freeNodes = p.level(end).enum.freeNodes;
%x0(freeNodes) = x;
%p.level(end).x0 = x0;
%p.level(end).x = x0;

DWRx0 = p.level(end).DWRx0;
freeNodes = p.level(end).enum.freeNodes;
DWRx0(freeNodes) = DWRx;
p.level(end).DWRx0 = DWRx0;
p.level(end).DWRx = DWRx0;

postProc = p.statics.postProc;
p = postProc(p);

pdeSolver = p.params.pdeSolver;
%% get function-value E(x)
getFuncVal = str2func(['TWP1DWRDWRgetFuncValPR']);
p = getFuncVal(p);
val = p.level(end).funcVal;

%% get the jacobian of E(x)
if nargout > 1
    getJacobian = str2func(['TWP1DWRDWRgetJacobianPR']);
    p = getJacobian(p);
    jacobian = p.level(end).jacobi;
end
