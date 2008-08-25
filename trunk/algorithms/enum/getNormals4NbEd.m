function normals4NbEd = getNormals4NbEd(c4n,Nb)
% normals4NbEd
% The function getNormals4NbEd returns a nrNbEdges times 2 matrix. The input 
% is c4n and Nb. The matrix contains in each row j the two coordinates of the 
% outer unit normal at the j-th Neumann edge, corresponding to the order in NbEdges.

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


normals4NbEd = [];
if(~isempty(Nb))
	tangents = ( c4n(Nb(:,2),:) - c4n(Nb(:,1),:) );
	dummy = ones(size(tangents,1),2);
	normals = tangents(:,[2 1]).*[dummy(:,1),-dummy(:,2)];
	lengthNormal = sqrt(normals(:,1).^2 + normals(:,2).^2);

	normals4NbEd  = normals./[lengthNormal, lengthNormal];
end
