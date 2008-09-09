function angles4e = getAngles4e(tangents4e)
% angles4e
% The function getAngles4e returns a [3 nrElems] matrix.
% The input is tangents4e produced by the function getTangente4e.
% The matrix contains in each column $j$ the three interior angles of element 
% j. The order of the angles correspond to the order of the nodes in n4e,
% i.e. the angles at the first node of each element are in the first row,
% the one at the second in the second etc.

% Copyright 2007 Jan Reininghaus, David Guenther
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


tangent4eU = tangents4e(:,1,:);
tangent4eV = tangents4e(:,2,:);

dummyU = -tangent4eU .* tangent4eU([3 1 2],:,:);
dummyV = -tangent4eV .* tangent4eV([3 1 2],:,:);

dummy = dummyU + dummyV;
dummy = reshape(dummy,3,[]);
angles4e = acos(dummy)';
