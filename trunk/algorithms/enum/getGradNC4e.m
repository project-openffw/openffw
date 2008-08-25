function gradNC4e = getGradNC4e(midpoint4ed,ed4e,area4e)
% gradNC4e
% The function getGrad4e returns a [3 2 nrElems] matrix. 
% The input is c4n n4e and area4e produced by the function getArea4e. 
% The three dimensional matrix contains in the j-th [3 2] matrix 
% the coordinates of the local gradients of the three P1-NC basis function
% which are not constantly zero on the j-th element. 
% The first row in each [3 2] matrix contains the local gradient of the basis 
% function corresponding to the first edge in ed4e, 
% the second row the one of the second edge etc.

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


xCoord4e = reshape(midpoint4ed(ed4e,1),[],3);
yCoord4e = reshape(midpoint4ed(ed4e,2),[],3);

gradNC4eY =  yCoord4e(:,[2,3,1]) - yCoord4e(:,[3,1,2]);
gradNC4eX =  xCoord4e(:,[3,1,2]) - xCoord4e(:,[2,3,1]);

gradNC4e = zeros(2*size(gradNC4eY,1),size(gradNC4eY,2));
gradNC4e(1:2:end,:) = gradNC4eY;
gradNC4e(2:2:end,:) = gradNC4eX;

dummy = repmat(area4e',2,1)/4;
dummy = repmat(dummy(:),1,3);
gradNC4e = gradNC4e./(2*dummy);

gradNC4e = reshape(gradNC4e',3,2,[]);
