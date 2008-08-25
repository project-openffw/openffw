function area4n = getArea4n(area4e,e4n)
% area4n
% The function getArea4n returns a nrNodes times 1 matrix. The input is
% e4n and area4e produced by the functions getE4n and getArea4e,
% respectively. The matrix contains in each row j the area of the
% patch of node j.

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


[I,J,S] = find(e4n);
newS = area4e(S);
area4n = sparse(I,J,newS);
area4n = full(sum(area4n,2));
