function triplot2(n4e,x,y)
% plots a triangular mesh in 2D
%
% Usage: triplot2(n4e,c4n(:,1),c4n(:,2));
% n4e - nodes for elements:    each triangle corresponds to one row of three vertices
% c4n - coordinates for nodes: each vertex corresponds to one row of coordinates
%       c4n(:,1) consists of the x-coordinates of the vertices
%       c4n(:,2) consists of the y-coordinates of the vertices

% Copyright 2007 Joscha Gedicke, Hella Rabus
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

clf;
 hold on; 
 for i = 1 : size(n4e,1)
     curElem = n4e(i,[1 2 3 1]);
     plot(x(curElem),y(curElem),'b');
 end
 hold off

%patch('Faces',n4e,'Vertices',[x,y],'facecolor','white','edgecolor','b');
%patch(x(n4e)',y(n4e)','c');
