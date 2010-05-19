function [val,jacobian] = ODP1getNonLinearSolution(x,p)
%author: David Guenther
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
