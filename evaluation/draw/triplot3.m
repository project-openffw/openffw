function triplot3(n4ed,x,y,z)
% plots a triangular mesh in 2D (red triangles)
%       and its corresponding solution in 3D (blue triangles)
%
% Use: triplot3(n4e(:,[1 2 3 1]),c4n(:,1),c4n(:,2),u)
%  or  triplot3(n4ed,c4n(:,1),c4n(:,2),u)
%  n4e - nodes for elements: each triangle corresponds to one row of three vertices.
%  c4n - coordinates for nodes: each vertex corresponds to one row of coordinates.
%        x=c4n(:,1) consists of the x-coordinates of the vertices
%        y=c4n(:,2) consists of the y-coordinates of the vertices
%   z=u  - is the solution, it consists of the discrete solution at each vertex
%
% n4ed - nodes for edges: each edge corresponds to one row of two vertices.

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

hold on; 
for i = 1 : size(n4ed,1)
    curEdge =n4ed(i,:);
    plot3(x(curEdge),y(curEdge),zeros(size(z(curEdge))),'r');
    plot3(x(curEdge),y(curEdge),z(curEdge),'b');
end
hold off


for d=-20:10:80
    view ( d, 30 );
    pause(1)
end
