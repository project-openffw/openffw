function p = ODP2enumerate(p)
% author: David Guenther 
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

%% get generic informations
p = genericEnumerate(p);

%% INPUT 
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Db = p.level(end).geom.Db;
area4e = p.level(end).enum.area4e;
ed4e = p.level(end).enum.ed4e;
DbEd = p.level(end).enum.DbEd;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;

%% get discretization-specific informations
fixedNodes = [unique(Db);DbEd+nrNodes];
freeNodes = setdiff(1:nrNodes+nrEdges, fixedNodes);

% Nr. of Degrees of Freedom
nrDoF = length(freeNodes);

grad4e = getGrad4e(c4n,n4e,area4e);

dofU4e = [n4e, ed4e+nrNodes];

%% OUTPUT
p.level(end).enum.dofU4e = dofU4e;
p.level(end).nrDoF = nrDoF;
p.level(end).enum.P1grad4e = grad4e;
p.level(end).enum.freeNodes = freeNodes;
p.level(end).enum.fixedNodes = fixedNodes;
