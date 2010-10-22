function [val,jacobian] = TWP1DWRgetNonLinearSolution(DWRx,p)
% author: David Guenther, Lena Noack

%% save the variable DWRx
%p.level(end).x0DWR = DWRx;
%p.level(end).xDWR = DWRx;

DWRx0 = p.level(end).x0DWR;
freeNodes = p.level(end).enum.freeNodes;
DWRx0(freeNodes) = DWRx;
p.level(end).x0DWR = DWRx0;
p.level(end).xDWR = DWRx0;
%p.level(end).xDWR = p.level(end).x;


postProc = p.statics.postProc;
p = postProc(p);

pdeSolver = p.params.pdeSolver;
%% get function-value E(x)
%getFuncVal = str2func(['TWP1getFuncVal']);
getFuncVal = str2func(['TWP1DWRgetFuncValIM']);
p = getFuncVal(p);
val = p.level(end).funcVal;

%% get the jacobian of E(x)
if nargout > 1
%    getJacobian = str2func(['TWP1getJacobian']);
    getJacobian = str2func(['TWP1DWRgetJacobianIM']);
    p = getJacobian(p);
    jacobian = p.level(end).jacobi;
end
