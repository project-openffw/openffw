function p = ODP1prolong(p,lvl)
%author: David Guenther
% Copyright 2007 David Guenther
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
%
%
%

%% INPUT
Db = p.level(lvl).geom.Db;
c4n = p.level(lvl).geom.c4n;
n4e = p.level(lvl).geom.n4e;
parents4e = p.level(lvl).enum.parents4e;
nrNodes = p.level(lvl).nrNodes;
nrElems = p.level(lvl).nrElems;

u_D = p.problem.u_D;
u_h = p.statics.u_h;

% x0 = p.level(1).P1x;

%% prolongation
DbNodes = unique(Db);
x0 = ones(nrNodes,1);

if lvl > 2
    x0 = zeros(nrNodes,1);
    for curElem = 1:nrElems
        nodes = n4e(curElem,:);
        coords = c4n(nodes,:);
        parent = parents4e(curElem);

        prolongU = u_h(coords,parent,lvl-1,p);
        x0(nodes) = prolongU;
    end
end

x0(unique(Db)) = u_D(c4n(DbNodes,:),p);

%% OUTPUT
p.level(lvl).x = x0;
