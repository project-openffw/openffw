function [val,jacobian] = ODRTGOALgetNonLinearSolution(x,p)
% author: David Guenther
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
