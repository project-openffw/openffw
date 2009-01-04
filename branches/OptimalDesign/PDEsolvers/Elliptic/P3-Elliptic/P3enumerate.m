function p = P3enumerate(p)
% creates all necessarily data 
% for a conforming P3-FE method.

% Copyright 2007 Joscha Gedicke
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


%% Generic Enumeration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = genericEnumerate(p);

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Db = p.level(end).geom.Db;
area4e = p.level(end).enum.area4e;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;
DbEd = p.level(end).enum.DbEd;
ed4e = p.level(end).enum.ed4e;

%% Enumeration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradients of P1-Hat functions
P1grad4e = getGrad4e(c4n,n4e,area4e);

% Nodes on Dirichlet boundary
fixedNodes = [unique(Db);nrNodes+DbEd;nrNodes+nrEdges+DbEd];
freeNodes = setdiff(1:(nrNodes+2*nrEdges+nrElems), fixedNodes);

% Nr. of Degrees of Freedom
nrDoF = length(freeNodes);

% Degrees of Freedom for elements
dummy = [ zeros(nrEdges,1) , ones(nrEdges,1)];
dof4Ed = zeros(nrElems,6);
for curElem = 1 : nrElems
    curEdges = ed4e(curElem,:);
    dof4Ed(curElem,:) = ...
      [nrNodes+dummy(curEdges,1)'*nrEdges+curEdges, nrNodes+dummy(curEdges,2)'*nrEdges+curEdges];
    dummy(curEdges,:) = mod(dummy(curEdges,:)+1,2);
end
dofU4e = [n4e,dof4Ed,(nrNodes+2*nrEdges)+[1:nrElems]'];

%% Enumeration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).enum.dofU4e = dofU4e;
p.level(end).enum.P1grad4e = P1grad4e;
p.level(end).nrDoF = nrDoF;
p.level(end).enum.freeNodes = freeNodes;
p.level(end).enum.fixedNodes = fixedNodes;