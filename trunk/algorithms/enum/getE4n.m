function e4n = getE4n(n4e)
% e4n
% The function getE4n returns a nrNodes times nrNodes sparse matrix. 
% NOTE:
% The sparsity constant is bounded due to the used mesh generation. The
% input is n4e. In this matrix each entry (j,k) is the number of the
% element, whose boundary contains the nodes j and k in
% counter clockwise order as vertices or zero if the there is no such
% element. Since the nodes in n4e are oriented counter clockwise e4n
% gives you the number of the row in which the sequence j k is found. 
% Note that in this context k i j also contains the sequence j k.
% To find the patch of a node k, i.e all elements containing node k,
% just get the non zero entries of the k-th row or column.

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


nrElems = size(n4e,1);
I = [n4e(:,1);n4e(:,2);n4e(:,3)];
J = [n4e(:,2);n4e(:,3);n4e(:,1)];
S = [1:nrElems,1:nrElems,1:nrElems];
e4n = sparse(I,J,S);
