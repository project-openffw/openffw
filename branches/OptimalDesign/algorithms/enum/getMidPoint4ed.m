function midpoint4ed = getMidPoint4ed(c4n,n4ed)
% midPoint4ed
% The function getMidpoint4ed returns a [nrEdges 2] matrix.
% The input is c4n and n4ed. The matrix contains in each row j 
% the coordinates of the midpoint of the j-th edge.

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


coordX = c4n(:,1);
coordY = c4n(:,2);

coordX4n4ed = coordX(n4ed);
coordY4n4ed = coordY(n4ed);

midpointX = sum(coordX4n4ed,2)/2;
midpointY = sum(coordY4n4ed,2)/2;

midpoint4ed = [midpointX,midpointY];
