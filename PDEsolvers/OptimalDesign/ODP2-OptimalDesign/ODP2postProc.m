function p = ODP2postProc(p)
% author: David Guenther 

%% INPUT 
% load computed solution
x = p.level(end).x;

% load enumerated data
dofU4e = p.level(end).enum.dofU4e;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;

%% post-processing of calculated data
u = x(1:nrNodes+nrEdges);
u4e = zeros(nrElems,6);

for curElem = 1:nrElems
	localDoF = dofU4e(curElem,:);    
    curU = u(localDoF);
    u4e(curElem,:) = curU;
end

%% OUTPUT 
p.level(end).u = u;
p.level(end).u4e = u4e;
