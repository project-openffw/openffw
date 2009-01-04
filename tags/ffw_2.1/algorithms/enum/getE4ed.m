function e4ed = getE4ed(n4ed,e4n)
% e4ed
% The function getE4ed returns a nrEdges times 2 matrix. The input  is e4n
% and n4ed produced by the functions getE4n and getN4ed, respectively.
% The matrix contains in each row $j$ the element numbers of the elements
% which share the edge. If the edge is a boundary edge the second
% entry is zero.

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


e4ed = rowaddr(e4n,n4ed,n4ed(:,[2 1]));
e4ed = sort(e4ed,2,'descend');
