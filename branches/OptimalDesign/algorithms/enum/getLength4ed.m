function length4ed = getLength4ed(c4n,n4ed)
% length4ed
% The function getLength4ed returns a nrEdges times 1 matrix. The input is
% c4n and n4ed produced by the function getN4ed. The matrix contains in
% each row j the length of the $j$-th edge, according to the edge numbers
% created in getEd4n.

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

coordX4edge = coordX(n4ed);
coordY4edge = coordY(n4ed);

xLength = coordX4edge(:,2) - coordX4edge(:,1);
yLength = coordY4edge(:,2) - coordY4edge(:,1);

length4ed = sqrt(xLength.^2 + yLength.^2);
