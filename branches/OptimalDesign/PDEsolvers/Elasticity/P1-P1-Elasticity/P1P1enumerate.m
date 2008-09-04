function p = P1P1enumerate(p)
%enumerate.m creates all necessarily data for a P1-FE method 
%in linear elasticity. 
%
%authors: David Guenther, Jan Reininghaus

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