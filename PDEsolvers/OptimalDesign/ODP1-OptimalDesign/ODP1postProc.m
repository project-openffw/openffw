function p = ODP1postProc(p)
%author: David Guenther, Lena Noack

%% INPUT 
% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;

% load computed solution
x = p.level(end).x;
DWRx = loadField('p.level(end)','DWRx',p,x);

% load enumerated data
P1grad4e = p.level(end).enum.grad4e;
nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
area4n = p.level(end).enum.area4n;
area4e = p.level(end).enum.area4e;

wh = p.statics.basisU;
qh = p.statics.stressBasis;

nonLinearExactDer = p.problem.nonLinearExactDer;
lvl = size(p.level,2);
%% post-processing of calculated data
u = x(1:nrNodes);
DWRu = DWRx(1:nrNodes);

u4e = zeros(nrElems,3);
DWRu4e = zeros(nrElems,3);
grad4e = zeros(nrElems,2);
DWRgrad4e = zeros(nrElems,2);
Aph = zeros(nrNodes,2);
AGradh = zeros(nrNodes,2);
DWRAGradh = zeros(nrNodes,2);
Auh = zeros(nrNodes,1);
DWRAuh = zeros(nrNodes,1);

Awh = zeros(nrNodes,1);
AGrad_wh = zeros(nrNodes,2);

for curElem = 1:nrElems
	curGrads = P1grad4e(:,:,curElem);
	curNodes = n4e(curElem,:);
	area = area4e(curElem);
    
    curU = u(curNodes);
    curDWRU = DWRu(curNodes);
    u4e(curElem,:) = curU;
    DWRu4e(curElem,:) = curDWRU;

    grad4e(curElem,:) = curU'*curGrads;
    DWRgrad4e(curElem,:) = curDWRU'*curGrads;
    curSigma = nonLinearExactDer(norm(grad4e(curElem,:)),curElem,lvl,p)*grad4e(curElem,:);
    Aph(curNodes,:) = Aph(curNodes,:) + (area*curSigma'*[1,1,1])';
    AGradh(curNodes,:) = AGradh(curNodes,:) + (area*grad4e(curElem,:)'*[1,1,1])';
    DWRAGradh(curNodes,:) = DWRAGradh(curNodes,:) + (area*DWRgrad4e(curElem,:)'*[1,1,1])';
    Auh(curNodes,:) = Auh(curNodes,:) + area*curU'*[1;1;1];
    DWRAuh(curNodes,:) = DWRAuh(curNodes,:) + area*curDWRU'*[1;1;1];

%    curWh = wh(c4n(curNodes,1),c4n(curNodes,2),curElem,lvl,p)'*[1;1;1]; %???? 3x1
    curWh = wh(c4n(curNodes,1),c4n(curNodes,2),curElem,lvl,p)'; %3xNodes 
    curQh = qh(c4n(curNodes,1),c4n(curNodes,2),curElem,lvl,p); %3x2xNodes
%    gradWh4e(curElem,:) = curWh'*curGrads;
    gradWh4e = matMul(reshape(curWh,[1 3 3]),curQh); %1x2xNodes
    Awh(curNodes,:) = Awh(curNodes,:) + area*curWh'*[1;1;1];
%    Awh(curNodes,:) = Awh(curNodes,:) + (area*[1 1 1]*curWh)';
%    AGrad_wh(curNodes,:) = AGrad_wh(curNodes,:) + (area*gradWh4e(curElem,:)'*[1,1,1])';
    AGrad_wh(curNodes,:) = AGrad_wh(curNodes,:) + reshape(area*gradWh4e,[3 2]);  %3x2
end

area4nVec = area4n*[1,1];
Aph = Aph./area4nVec;
AGradh = AGradh./area4nVec;
DWRAGradh = DWRAGradh./area4nVec;
Auh = Auh./area4n;
DWRAuh = DWRAuh./area4n;

AGrad_wh = AGrad_wh./area4nVec;
Awh = Awh./area4n;

%% OUTPUT 
p.level(end).u = u;
p.level(end).DWRu = DWRu;
p.level(end).u4e = u4e;
p.level(end).DWRu4e = DWRu4e;
p.level(end).grad4e = grad4e;
p.level(end).DWRgrad4e = DWRgrad4e;
p.level(end).Aph = Aph;
p.level(end).AGradh = AGradh;
p.level(end).DWRAGradh = DWRAGradh;
p.level(end).Auh = Auh;
p.level(end).DWRAuh = DWRAuh;

p.level(end).AGrad_wh = AGrad_wh;
p.level(end).Awh = Awh;
