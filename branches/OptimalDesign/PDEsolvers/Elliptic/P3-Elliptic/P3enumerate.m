function p = P3enumerate(p)

% author: Joscha Gedicke

p = genericEnumerate(p);

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Db = p.level(end).geom.Db;
area4e = p.level(end).enum.area4e;
nrNodes = p.level(end).nrNodes;
n4ed    = p.level(end).enum.n4ed;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;
DbEd = p.level(end).enum.DbEd;
ed4e = p.level(end).enum.ed4e;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).enum.dofU4e = dofU4e;
p.level(end).enum.P1grad4e = P1grad4e;
p.level(end).nrDoF = nrDoF;
p.level(end).enum.freeNodes = freeNodes;
p.level(end).enum.fixedNodes = fixedNodes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
