function [val,jacobian] = TWDWRgetNonLinearSolution(x,p)
% author: Lena Noack

%% save the variable x
p.level(end).x0 = x;
p.level(end).x = x;

p = TWDWRpostProc(p);

%% get function-value E(x)
p = TWDWRgetFuncVal(p);
val = p.level(end).funcVal;
p = TWDWRgetJacobian(p);
%% get the jacobian of E(x)
if nargout > 1
    p = TWDWRgetJacobian(p);
    jacobian = p.level(end).jacobi;
end
