function p = TWDWRpostProc(p)
% author: Lena Noack
%postproc.m computes the post-processing datas for 
%the nonlinear mixed RT0-P0-FE method.

%% INPUT 
% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;

% load computed solution
x = p.level(end).x;

% load enumerated data
ed4e = p.level(end).enum.ed4e;
area4e = p.level(end).enum.area4e;
area4n = p.level(end).enum.area4n;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;
nrNodes = p.level(end).nrNodes;

stressBasis = p.statics.stressBasis;

lvl = size(p.level,2);

%% post-processing of calculated data
u = x(nrEdges+1:nrEdges+nrElems)';
gradP = x(1:nrEdges);
lambda_1 = x(nrEdges+nrElems+1 : nrEdges+nrElems+nrEdges);
lambda_2 = x(nrEdges+nrElems+nrEdges+1:end);

grad = zeros(3,2,nrElems);
lambda1 = zeros(3,2,nrElems);
grad4e = zeros(nrElems,3);
lambda14e = zeros(nrElems,3);

Aph = zeros(nrNodes,2);
Alambda1 = zeros(nrNodes,2);
Auh = zeros(nrNodes,1);
Alambda2 = zeros(nrNodes,1);

for curElem = 1:nrElems
    curEdges = ed4e(curElem,:);
	curNodes = n4e(curElem,:);
	curCoords = c4n(curNodes,:)';
	area = area4e(curElem);
    
    curP = x(curEdges);
    curLambda1 = lambda_1(curEdges);
    
    grad4e(curElem,:) = curP;
    lambda14e(curElem,:) = curLambda1;
    
    P1 = curCoords(:,1);
    P2 = curCoords(:,2);
    P3 = curCoords(:,3);
    
    evalBasisP1 = stressBasis(P1(1),P1(2),curElem,lvl,p);
    evalBasisP2 = stressBasis(P2(1),P2(2),curElem,lvl,p);
    evalBasisP3 = stressBasis(P3(1),P3(2),curElem,lvl,p);
    
    grad(1,:,curElem) = curP'*evalBasisP1;
    grad(2,:,curElem) = curP'*evalBasisP2;
    grad(3,:,curElem) = curP'*evalBasisP3;
    
    lambda1(1,:,curElem) = curLambda1'*evalBasisP1;
    lambda1(2,:,curElem) = curLambda1'*evalBasisP2;
    lambda1(3,:,curElem) = curLambda1'*evalBasisP3;
    
    Aph(curNodes,:) = Aph(curNodes,:) + area*grad(:,:,curElem);
    Alambda1(curNodes,:) = Alambda1(curNodes,:) + area*lambda1(:,:,curElem);
    
    Auh(curNodes) = Auh(curNodes) + area*u(curElem);
    Alambda2(curNodes) = Alambda2(curNodes) + area*lambda_2(curElem);
end

Auh = Auh./area4n;
Alambda2 = Alambda2./area4n;

area4n = area4n*[1,1];
Aph = Aph./area4n;
Alambda1 = Alambda1./area4n;

%% OUTPUT 
p.level(end).u = u;
p.level(end).u4e = u;
p.level(end).Auh = Auh;

p.level(end).lambda_2 = lambda_2;
p.level(end).lambda24e = lambda_2;
p.level(end).Alambda2 = Alambda2;

p.level(end).gradP = gradP;
p.level(end).grad = grad;
p.level(end).grad4e = grad4e;
p.level(end).Aph = Aph;

p.level(end).lambda_1 = lambda_1;
p.level(end).lambda1 = lambda1;
p.level(end).lambda14e = lambda14e;
p.level(end).Alambda1 = Alambda1;
