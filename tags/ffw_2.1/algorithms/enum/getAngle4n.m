function angle4n = getAngle4n(angles4e,n4e,nrElems,nrNodes)
% angle4n
% The function getAngle4n returns a [nrNodes 1] matrix. 
% The input is angles4e, n4e, nrElems and nrNodes where angles4e is generated
% by getAngles4e. The matrix contains in each row j the angle at the node j,
% where the angle is 2 pi for inner nodes and the angle inside the domain 
% and between the two boundary edges containing the node for nodes at the boundary.

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


angle4n = zeros(nrNodes,1);

for curElem = 1:nrElems
	curNodes = n4e(curElem,:);
	curAngles = angles4e(curElem,:);
	angle4n(curNodes) = angle4n(curNodes) + curAngles';
end
	
