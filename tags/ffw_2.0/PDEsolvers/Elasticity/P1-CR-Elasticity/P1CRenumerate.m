function p = P1CRenumerate(p)
% creates all necessarily data 
% for the Kouhia-Stenberg FE in linear elasticity.

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


p = genericEnumerate(p);

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Db = p.level(end).geom.Db;
DbEd = p.level(end).enum.DbEd;
midPoint4ed = p.level(end).enum.midPoint4ed;
ed4e = p.level(end).enum.ed4e;
area4e = p.level(end).enum.area4e;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gradients of P1-Hat functions
grad4e = getGrad4e(c4n,n4e,area4e);

% gradients of P1NC-Hat functions
gradNC4e = getGradNC4e(midPoint4ed,ed4e,area4e);

fixedNodes = unique(Db);
freeNodes = setdiff(1:nrNodes, fixedNodes);

fixedEdges = DbEd;
freeEdges = setdiff(1:nrEdges, fixedEdges);

nrDoF = length(freeNodes) + length(freeEdges);

dofU4e =  [n4e,nrNodes + ed4e];
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).enum.dofU4e = dofU4e;
p.level(end).enum.grad4e = grad4e;
p.level(end).enum.gradNC4e = gradNC4e;
p.level(end).enum.freeNodes = freeNodes;
p.level(end).enum.fixedNodes = fixedNodes;
p.level(end).enum.freeEdges = freeEdges;
p.level(end).enum.fixedEdges = fixedEdges;
p.level(end).nrDoF = nrDoF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
