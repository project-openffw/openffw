function p = ODRTpostProc(p)
% author: David Guenther, Lena Noack
%postproc.m computes the post-processing datas for
%the nonlinear mixed RT0-P0-FE method.

%% INPUT
% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;

% load computed solution
x = p.level(end).x;
DWRx = loadField('p.level(end)','DWRx',p,x);
%DWRx = p.level(end).DWRx;

% load enumerated data
P1grad4e = p.level(end).enum.grad4e;
ed4e = p.level(end).enum.ed4e;
area4e = p.level(end).enum.area4e;
area4n = p.level(end).enum.area4n;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;
nrNodes = p.level(end).nrNodes;

stressBasis = p.statics.stressBasis;

lvl = size(p.level,2);

%% post-processing of calculated data
u = x(nrEdges+1:end)';          %solution u
DWRu = DWRx(nrEdges+1:end)';    %dual solution z

gradU4e = zeros(nrElems,2);   %primal gradient Du
DWRgradU4e = zeros(nrElems,2);%dual gradient Du

grad = zeros(3,2,nrElems);
DWRgrad = zeros(3,2,nrElems);
grad4e = zeros(nrElems,3);    %primal sigma
DWRgrad4e = zeros(nrElems,3); %dual sigma

Aph = zeros(nrNodes,2);
Agradh = zeros(nrNodes,2);
DWRAph = zeros(nrNodes,2);
Auh = zeros(nrNodes,1);
DWRAuh = zeros(nrNodes,1);

for curElem = 1:nrElems
    curEdges = ed4e(curElem,:);
    curNodes = n4e(curElem,:);
    curCoords = c4n(curNodes,:)';
    area = area4e(curElem);

	curGrads = P1grad4e(:,:,curElem);
    curU = u(curNodes);
    u4e(curElem,:) = curU;
    curDWRU = DWRu(curNodes);
    DWRu4e(curElem,:) = curDWRU;
    DWRgradU4e(curElem,:) = (curDWRU*curGrads);
    gradU4e(curElem,:) = curU*curGrads;

    curP = x(curEdges);             %sigma
    DWRcurP = DWRx(curEdges);       %sigma (dual solution)
    grad4e(curElem,:) = curP;
    DWRgrad4e(curElem,:) = DWRcurP;

    P1 = curCoords(:,1);
    P2 = curCoords(:,2);
    P3 = curCoords(:,3);

    evalBasisP1 = stressBasis(P1(1),P1(2),curElem,lvl,p);
    evalBasisP2 = stressBasis(P2(1),P2(2),curElem,lvl,p);
    evalBasisP3 = stressBasis(P3(1),P3(2),curElem,lvl,p);

    grad(1,:,curElem) = curP'*evalBasisP1;
    grad(2,:,curElem) = curP'*evalBasisP2;
    grad(3,:,curElem) = curP'*evalBasisP3;

    DWRgrad(1,:,curElem) = DWRcurP'*evalBasisP1;
    DWRgrad(2,:,curElem) = DWRcurP'*evalBasisP2;
    DWRgrad(3,:,curElem) = DWRcurP'*evalBasisP3;

    Aph(curNodes,:) = Aph(curNodes,:) + area*grad(:,:,curElem);
    Auh(curNodes) = Auh(curNodes) + area*u(curElem);
    DWRAph(curNodes,:) = DWRAph(curNodes,:) + area*DWRgrad(:,:,curElem);
    DWRAuh(curNodes) = DWRAuh(curNodes) + area*DWRu(curElem);
    Agradh(curNodes,:) = Agradh(curNodes,:) + (area*gradU4e(curElem,:)'*[1,1,1])';
end

Auh = Auh./area4n;
DWRAuh = DWRAuh./area4n;
area4n = area4n*[1,1];
Aph = Aph./area4n;
Agradh = Agradh./area4n;
DWRAph = DWRAph./area4n;

%% OUTPUT
p.level(end).u = u;
p.level(end).u4e = u4e;%u;
p.level(end).gradU4e = gradU4e;
p.level(end).grad4e = grad4e;
p.level(end).Aph = Aph;
p.level(end).Agradh = Agradh;
p.level(end).Auh = Auh;

p.level(end).DWRu = DWRu;
p.level(end).DWRu4e = DWRu;
p.level(end).DWRgrad4e = DWRgrad4e;
p.level(end).DWRgradU4e = DWRgradU4e;
p.level(end).DWRAph = DWRAph;
p.level(end).DWRAuh = DWRAuh;

%% get projection error
% q = p;
% q = ODRTestimate_Jump(q);
% p.level(end).errorJump = q.level(end).estimatedError;

%% compute discrete conjugate energy
% conjEnergy = sum(integrate(n4e,lvl,10,@getConjEnergy,p));
% p.level(end).conjEnergy = conjEnergy;

%% compute area contact zone
% volumeFraction = p.statics.volumeFraction;
%  
% coordsX = reshape(c4n(n4e,1),[],3)';
% coordsY = reshape(c4n(n4e,2),[],3)';
% 
% val = zeros(3,nrElems);
% for curElem = 1:nrElems
%     val(:,curElem) = volumeFraction(coordsX(:,curElem),coordsY(:,curElem),curElem,lvl,p);
% end
% 
% dummy = sum(val,1)';
% I = find((dummy ~= 3) & (dummy ~= 0));
% totalArea = sum(area4e(I));
% p.level(end).areaContactZone = totalArea;

%% compute average error
% averageError = integrate(n4e,lvl,10,@phminusAph,p);
% p.level(end).phminusAph = sqrt(sum(averageError));
% p.level(end).phminusAph4e = sqrt(averageError);

%% supply the error ||p_{h,\epsilon)-Ap_{h,\epsilon)||_L^2
function val = phminusAph(x,y,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
Asigma_h = p.statics.Ap_h;

Asigmah = Asigma_h(x,y,curElem,lvl,p);
sigmah = sigma_h(x,y,curElem,lvl,p);

val = sum( (sigmah - Asigmah).*(sigmah - Asigmah),2 );
val = reshape(val,[1 1 length(x)]);

%% supply conjugate energy functional
function evalEnergy = getConjEnergy(x,y,curElem,lvl,p)

conjFunctional = p.problem.nonLinearConjExact;
sigma_h = p.statics.sigma_h;

evalSigma = sigma_h(x,y,curElem,lvl,p);
absSigma = (evalSigma(:,1).^2 + evalSigma(:,2).^2).^(1/2);

evalEnergy = conjFunctional(absSigma,curElem,lvl,p);

evalEnergy = reshape(evalEnergy,[1 1 length(x)]);
