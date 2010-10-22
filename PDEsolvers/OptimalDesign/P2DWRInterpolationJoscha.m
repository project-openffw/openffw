function p = P2DWRInterpolation(p,curLvl)

% Copyright 2007 Joscha Gedicke, Lena Noack
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

if nargin < 2
    curLvl = length(p.level);
end

%% Input

freeNodes = p.level(end).enum.freeNodes;
midPoint4ed = p.level(end).enum.midPoint4ed;
ed4e = p.level(end).enum.ed4e;
e4n = p.level(end).enum.e4n;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
c4n = p.level(end).geom.c4n;
n4e = p.level(end).geom.n4e;
uh = p.level(end).DWRu;
n4ed = p.level(end).enum.n4ed;

%% P2 Interpolation

u = zeros(nrNodes+nrEdges,1);

[I,J,S] = find(e4n');

for curNode = freeNodes
    [a,b] = binsearch(J,curNode);
    patch = S(a:b);
    edges = ed4e(patch,:);
%     edges = unique(ed4e(patch,:));
    x = midPoint4ed(edges,1);
    y = midPoint4ed(edges,2);
    A = [ ones(length(x),1), x, y, x.*y, x.^2, y.^2];
    b = 1/2*(uh(n4ed(edges,1))+uh(n4ed(edges,2)));
    c = A\b;
    x = c4n(curNode,1);
    y = c4n(curNode,2);
    u(curNode) = [1,x, y, x.*y, x.^2, y.^2]*c;
end

u(nrNodes+1:end) = 1/2*( u(n4ed(:,1)) + u(n4ed(:,2)) );

dofU4e = [n4e,nrNodes+ed4e];

%% Output
p.level(curLvl).P2DWRu4e = u(dofU4e);
