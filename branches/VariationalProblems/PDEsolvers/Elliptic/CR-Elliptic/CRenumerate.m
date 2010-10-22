function p = CRenumerate(p)
%enumerate.m creates all necessarily data for a 
%nonconforming CR-FE method.
%
%authors: David Guenther, Jan Reininghaus
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
