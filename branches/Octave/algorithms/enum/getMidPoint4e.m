function midPoint4e = getMidPoint4e(c4n,n4e)
% midPoint4e
% The function getMidpoint4e returns a [nrElems 2] matrix. 
% The input is c4n and n4e. The matrix contains in each row j
% the coordinates of the midpoint of element j.

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


coord4n4e = reshape( c4n(n4e',:),3,[] );
s = sum(coord4n4e,1);
midPoint4e = reshape(s,[],2)/3;
