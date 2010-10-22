function p = ODRTGOALenumerate(p)
%enumerate.m creates all necessarily data for the nonlinear - mixed
%RT0-P0-FE method.
% author: Lena Noack
%% get generic informations
p = genericEnumerate(p);

%% INPUT 
NbEd = p.level(end).enum.NbEd;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
ed4e = p.level(end).enum.ed4e;
area4e = p.level(end).enum.area4e;

%% get discretization-specific informations

fixedNodes = [NbEd,nrEdges+nrElems+NbEd];
freeNodes = setdiff(1:(2*nrEdges+2*nrElems), fixedNodes);

% Nr. of Degrees of Freedom
nrDoF = length(freeNodes);
grad4e = getGrad4e(c4n,n4e,area4e);

dofU4e = (1:nrElems)';
dofSigma4e = ed4e;

%% OUTPUT 
p.level(end).enum.dofU4e = dofU4e;
p.level(end).enum.dofSigma4e = dofSigma4e;
p.level(end).nrDoF = nrDoF;
p.level(end).enum.grad4e = grad4e;
p.level(end).enum.freeNodes = freeNodes;
p.level(end).enum.fixedNodes = fixedNodes;
