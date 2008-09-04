function [val,jacobian] = ODRTGOALgetNonLinearSolution(x,p)
% author: David Guenther

%% save the variable x
p.level(end).x0 = x;
p.level(end).x = x;

p = ODRTGOALpostProc(p);

%% get function-value E(x)
p = ODRTGOALgetFuncVal(p);
val = p.level(end).funcVal;
p = ODRTGOALgetJacobian(p);
%% get the jacobian of E(x)
if nargout > 1
    p = ODRTGOALgetJacobian(p);
    jacobian = p.level(end).jacobi;
end
