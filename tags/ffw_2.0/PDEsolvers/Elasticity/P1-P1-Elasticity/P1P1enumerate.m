function p = P1P1enumerate(p)
% creates all necessarily data for a P1-FE method 
% in linear elasticity. 

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
nrNodes = p.level(end).nrNodes;
area4e = p.level(end).enum.area4e;
Db = p.level(end).geom.Db;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gradients of P1-Hat functions
grad4e = getGrad4e(c4n,n4e,area4e);

fixedNodes = unique(Db);
freeNodes = setdiff(1:nrNodes, fixedNodes);

nrDoF = 2 * length(freeNodes);

dofU4e = [n4e,nrNodes + n4e];
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).enum.dofU4e = dofU4e;
p.level(end).enum.grad4e = grad4e;
p.level(end).enum.freeNodes = freeNodes;
p.level(end).enum.fixedNodes = fixedNodes;
p.level(end).nrDoF = nrDoF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
