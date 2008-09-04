function p = P1enumerate(p)
%enumerate.m creates all necessarily data for a 
%conforming P1-FE method.
%
%authors: David Guenther, Jan Reininghaus

p = genericEnumerate(p);

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Db = p.level(end).geom.Db;
area4e = p.level(end).enum.area4e;
nrNodes = p.level(end).nrNodes;
n4ed    = p.level(end).enum.n4ed;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gradients of P1-Hat functions
grad4e = getGrad4e(c4n,n4e,area4e);

% Nodes on Dirichlet boundary
fixedNodes = unique(Db);
freeNodes = setdiff(1:nrNodes, fixedNodes);

% Nr. of Degrees of Freedom
nrDoF = length(freeNodes);

dofU4e = n4e;

dofU4ed = n4ed;

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).enum.dofU4e = dofU4e;
p.level(end).enum.dofU4ed = dofU4ed;
p.level(end).enum.grad4e = grad4e;
p.level(end).nrDoF = nrDoF;
p.level(end).enum.freeNodes = freeNodes;
p.level(end).enum.fixedNodes = fixedNodes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
