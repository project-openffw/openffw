function p = P2enumerate(p)

% author: Joscha Gedicke

p = genericEnumerate(p);

%% INPUT
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Db = p.level(end).geom.Db;
area4e = p.level(end).enum.area4e;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
DbEd = p.level(end).enum.DbEd;
ed4e = p.level(end).enum.ed4e;

%% gradients of P1-Hat functions
P1grad4e = getGrad4e(c4n,n4e,area4e);

%% Nodes on Dirichlet boundary
fixedNodes = [unique(Db);nrNodes+DbEd];
freeNodes = setdiff(1:(nrNodes+nrEdges), fixedNodes);

%% Nr. of Degrees of Freedom
nrDoF = length(freeNodes);

%% Degrees of Freedom for elements
dofU4e = [n4e,nrNodes+ed4e];

%% OUTPUT
p.level(end).enum.dofU4e = dofU4e;
p.level(end).enum.P1grad4e = P1grad4e;
p.level(end).nrDoF = nrDoF;
p.level(end).enum.freeNodes = freeNodes;
p.level(end).enum.fixedNodes = fixedNodes;
