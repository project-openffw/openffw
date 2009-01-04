function p = ODP1postProc(p)
%author: David Guenther

%% INPUT 
% load geometry
n4e = p.level(end).geom.n4e;

% load computed solution
x = p.level(end).x;

% load enumerated data
P1grad4e = p.level(end).enum.grad4e;
nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
area4n = p.level(end).enum.area4n;
area4e = p.level(end).enum.area4e;

nonLinearExactDer = p.problem.nonLinearExactDer;
lvl = size(p.level,2);
%% post-processing of calculated data
u = x(1:nrNodes);

u4e = zeros(nrElems,3);
grad4e = zeros(nrElems,2);
Aph = zeros(nrNodes,2);

for curElem = 1:nrElems
	curGrads = P1grad4e(:,:,curElem);
	curNodes = n4e(curElem,:);
	area = area4e(curElem);
    
    curU = u(curNodes);
    u4e(curElem,:) = curU;
	grad4e(curElem,:) = curU'*curGrads;
    curSigma = nonLinearExactDer(norm(grad4e(curElem,:)),curElem,lvl,p)*grad4e(curElem,:);
    Aph(curNodes,:) = Aph(curNodes,:) + (area*curSigma'*[1,1,1])';
end

area4n = area4n*[1,1];
Aph = Aph./area4n;

%% OUTPUT 
p.level(end).u = u;
p.level(end).u4e = u4e;
p.level(end).grad4e = grad4e;
p.level(end).Aph = Aph;
