function area4e = getArea4e(c4n,n4e)
% area4e
% The function getArea4e returns a nrElems times 1 matrix. The input is
% n4e and c4n. The matrix contains in each row j the area of the
% element corresponding to the element number j.

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

area4e = 1/2 * ( (xCoord4e(:,2)-xCoord4e(:,1)).*(yCoord4e(:,3)-yCoord4e(:,1))...
		- (yCoord4e(:,2)-yCoord4e(:,1)).*(xCoord4e(:,3)-xCoord4e(:,1)) );
