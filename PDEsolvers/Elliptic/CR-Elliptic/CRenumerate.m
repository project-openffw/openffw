function p = CRenumerate(p)
% creates all necessarily data for a 
% nonconforming CR-FE method.

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
Db = p.level(end).geom.Db;
area4e = p.level(end).enum.area4e;
DbEd = p.level(end).enum.DbEd;
midPoint4ed = p.level(end).enum.midPoint4ed;
ed4e = p.level(end).enum.ed4e;
nrEdges = p.level(end).nrEdges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gradients of P1NC-Hat functions
gradNC4e = getGradNC4e(midPoint4ed,ed4e,area4e);

% Nodes on Dirichlet boundary
fixedNodes = unique(Db);

freeNodes = setdiff(1:nrEdges, DbEd);

% Nr. of Degrees of Freedom
nrDoF = length(freeNodes);

dofU4e = ed4e;

ed4ed = [1:nrEdges]';
dofU4ed = ed4ed;

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).enum.dofU4e = dofU4e;
p.level(end).enum.dofU4ed = dofU4ed;
p.level(end).enum.gradNC4e = gradNC4e;
p.level(end).nrDoF = nrDoF;
p.level(end).enum.freeNodes = freeNodes;
p.level(end).enum.fixedNodes = fixedNodes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
