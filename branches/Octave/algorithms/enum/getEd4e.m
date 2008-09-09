function ed4e = getEd4e(n4e, ed4n)
% ed4e
% The function getEd4e returns a nrElems \times 3 matrix. The input is n4e
% and ed4n produced by the function getEd4n. This matrix contains in
% each row j the number of the three edges of element j. The edge
% numbers are the ones generated  in the function GetEd4n. The edges
% are ordered counter clockwise beginning with the edge between the
% first and the second node in n4e.

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


ed4e = rowaddr(ed4n,n4e(:,[2 3 1]),n4e);
