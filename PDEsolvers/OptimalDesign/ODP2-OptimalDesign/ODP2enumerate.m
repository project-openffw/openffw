function p = ODP2enumerate(p)
% author: David Guenther 

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
