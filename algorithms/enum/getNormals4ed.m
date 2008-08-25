function normals4ed = getNormals4ed(normals4e,ed4e)
% normals4ed

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


normals4e = permute(normals4e, [3 2 1 ] );
normals4ed(ed4e(:,1),:) = normals4e(:,:,1); 
normals4ed(ed4e(:,2),:) = normals4e(:,:,2);
normals4ed(ed4e(:,3),:) = normals4e(:,:,3);
