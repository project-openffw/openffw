function grad4e = getGrad4e(c4n,n4e,area4e)
% grad4e
% The function getGrad4e returns a [3 2 nrElems] matrix. 
% The input is c4n n4e and area4e produced by the function getArea4e.
% The three dimensional matrix contains in the j-th [3 2] matrix 
% the coordinates of the local gradients of the three P1 basis function
% which are not constantly zero on the j-th element. 
% The first row in each [3 2] matrix contains the local gradient of the nodal basis 
% function for the first node of the element according to the order in n4e,
% the second row the one of the nodal basis function for the second node etc.

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


xCoord4e = reshape(c4n(n4e,1),[],3);
yCoord4e = reshape(c4n(n4e,2),[],3);

grad4eY =  yCoord4e(:,[2,3,1]) - yCoord4e(:,[3,1,2]);
grad4eX =  xCoord4e(:,[3,1,2]) - xCoord4e(:,[2,3,1]);

grad4e = zeros(2*size(grad4eY,1),size(grad4eY,2));
grad4e(1:2:end,:) = grad4eY;
grad4e(2:2:end,:) = grad4eX;

dummy = repmat(area4e',2,1);
dummy = repmat(dummy(:),1,3);
grad4e = grad4e./(2*dummy);

grad4e = reshape(grad4e',3,2,[]);
