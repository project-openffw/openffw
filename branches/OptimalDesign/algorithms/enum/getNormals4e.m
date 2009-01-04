function normals4e = getNormals4e(c4n,n4e,length4ed,ed4e)
% normals4e
% The function getNormnals4e returns a [3 2 nrElems] matrix. the input is 
% c4n, n4e and length4ed and ed4e produced by the functions getLength4ed
% and getEd4e. The three dimensional matrix contains in the j-th 3 times 2
% matrix the coordinates of the unit outer normals for the three edges of
% element j. The order of the three rows corresponds two the order in ed4e.

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

tangent4eX =  xCoord4e(:,[2,3,1]) - xCoord4e;
tangent4eY =  yCoord4e(:,[2,3,1]) - yCoord4e;

% Normalize 
lengthEd4e = reshape(length4ed(ed4e),[],3);
tangent4eX = tangent4eX ./ lengthEd4e;
tangent4eY = tangent4eY ./ lengthEd4e;

normals4e = zeros(2*size(tangent4eY,1),size(tangent4eY,2));
normals4e(1:2:end,:) = tangent4eY;
normals4e(2:2:end,:) = -tangent4eX;

normals4e = reshape(normals4e',3,2,[]);
