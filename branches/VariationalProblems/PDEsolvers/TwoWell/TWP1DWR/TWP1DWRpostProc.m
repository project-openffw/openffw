function p = P1postProc(p)
%author: David Guenther, Lena Noack

%% INPUT 
F1 = p.problem.F1;
F2 = p.problem.F2;

% load geometry
n4e = p.level(end).geom.n4e;

% load computed solution
x = p.level(end).x;
DWRx = loadField('p.level(end)','DWRx',p,x);

% load enumerated data
P1grad4e = p.level(end).enum.grad4e;
nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
area4n = p.level(end).enum.area4n;
area4e = p.level(end).enum.area4e;

CONV = p.params.CONV;
if strcmp(CONV,'c')
nonLinearExactDer1 = p.problem.nonLinearRegDerA;
nonLinearExactDer2 = p.problem.nonLinearRegDerB;
else
nonLinearExactDer1 = p.problem.nonLinearExactDerA;
nonLinearExactDer2 = p.problem.nonLinearExactDerB;
end

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

    evalNonLinearDer1 = nonLinearExactDer1(grad4e(curElem,1),grad4e(curElem,2),curElem,lvl,p);
    evalNonLinearDer2 = nonLinearExactDer2(grad4e(curElem,1),grad4e(curElem,2),curElem,lvl,p);
    if strcmp(CONV,'c')
    %% sigma = DW**(F) = W'**1(F)F + W'**2(F)(F2,F)F2
    curSigma = evalNonLinearDer1*grad4e(curElem,:) + ...
               evalNonLinearDer2*(grad4e(curElem,1)*F2(1) + grad4e(curElem,2)*F2(2))*F2';
    else
    %% sigma = DW(F) = W'1(F)(F-F1) + W'2(F)(F-F2)
    curSigma = evalNonLinearDer1*(grad4e(curElem,:)-F1') + evalNonLinearDer2*(grad4e(curElem,:)-F2');
    end

    Aph(curNodes,:) = Aph(curNodes,:) + (area*curSigma'*[1,1,1])';
    AGradh(curNodes,:) = AGradh(curNodes,:) + (area*grad4e(curElem,:)'*[1,1,1])';
    DWRAGradh(curNodes,:) = DWRAGradh(curNodes,:) + (area*DWRgrad4e(curElem,:)'*[1,1,1])';
    Auh(curNodes,:) = Auh(curNodes,:) + area*curU'*[1;1;1];
end

area4nVec = area4n*[1,1];
Aph = Aph./area4nVec;
AGradh = AGradh./area4nVec;
DWRAGradh = DWRAGradh./area4nVec;
Auh = Auh./area4n;

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
