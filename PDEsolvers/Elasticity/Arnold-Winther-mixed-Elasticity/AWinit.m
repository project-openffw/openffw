function p = AWinit(p)
% makes available all necessary initial data,
% handels the discrete displacement,stress, ...

% Copyright 2007 Jan Reininghaus, David Guenther
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.statics.basisU = @getDisplacementBasis;
p.statics.basisSigma = @getStressBasis;
p.statics.u_h = @getU_h;
p.statics.sigma_h = @getSigma_h;

% integration parameters 
%  -> up to which polynomial degree shall integration be exact?
p.params.integrationDegrees.createLinSys.Rhs = 9;


%% Basis Functions etc. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_h = getU_h(pts,curElem,lvl,p)
nrPts = size(pts,1);

u = p.level(lvl).u4e;
basisU = p.statics.basisU;
evalBasisU = basisU(pts,curElem,lvl,p);

nrBasisFuncU = size(p.level(end).enum.dofU4e,2);
curU = u(:,:,curElem);

curU = curU(:)';
curU = repmat(curU,2*nrPts,1);
curU = reshape(curU,[nrPts,2,nrBasisFuncU]);

u_h = evalBasisU.*curU;
u_hU = u_h(:,:,1:nrBasisFuncU/2);
u_hV = u_h(:,:,nrBasisFuncU/2+1:end);

u_hU = squeeze(sum(u_hU,3));
u_hV = squeeze(sum(u_hV,3));

u_h = [u_hU(:,1)'; u_hV(:,2)'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigma_h = getSigma_h(pts,curElem,lvl,p)
nrPts = size(pts,1);
sigma = p.level(lvl).sigma;

curSigma = sigma(curElem,:);
basisP3 = getBasisP3(pts,curElem,lvl,p);
basisCoefficents = p.level(lvl).basisCoefficents;
curBasisCoefficents = basisCoefficents(:,:,curElem);

nrPoints = nrPts/length(curElem);
curBasisCoefficents = permute(curBasisCoefficents,[2 1 3]);
curSigma = curSigma';
curSigma = reshape(curSigma,[24,1,length(curElem)]);
curSigma = permute(curSigma,[2 1 3]);
coeffs4sigma = matMul(curSigma,curBasisCoefficents);

I = (1:3:28);

% coeff a^(j)
coeffA = coeffs4sigma(:,I,:);
% coeff b^(j)
coeffB = coeffs4sigma(:,I+1,:);
% coeff c^(j)
coeffC = coeffs4sigma(:,I+2,:);

sigma_h1 = matMul(coeffA,basisP3);
sigma_h2 = matMul(coeffB,basisP3);
sigma_h3 = matMul(coeffC,basisP3);

%dimension sigma_h: 2x2 x nrEvaluationPoints x nrElems
sigma_h = zeros(2,2,nrPoints,length(curElem));
sigma_h(1,1,:,:) = sigma_h1;
sigma_h(2,2,:,:) = sigma_h2;
sigma_h(1,2,:,:) = sigma_h3;
sigma_h(2,1,:,:) = sigma_h3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function basisSigma = getStressBasis(x,y,curElem,lvl,p)
% 
% basisP3 = getBasisP3(x,y,curElem,lvl,p);
% basisCoefficents = p.level(lvl).basisCoefficents;
% curBasisCoefficents = basisCoefficents(:,:,curElem);
% 
% I = (1:3:28);
% 
% % coeff a^(j)
% coeffA = curBasisCoefficents(I);
% % coeff b^(j)
% coeffB = curBasisCoefficents(I+1);
% % coeff c^(j)
% coeffC = curBasisCoefficents(I+2);
% 
% basisSigma = [coeffA*basisP3; coeffC*basisP3;...
% 			coeffC*basisP3; coeffB*basisP3];
% 
% basisSigma = reshape(basisSigma,2,2,[]);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basisU = getDisplacementBasis(pts,curElem,lvl,p)
nrPts = size(pts,1);
x = pts(:,1);
y = pts(:,2);

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
area4e = p.level(lvl).enum.area4e;

n4e = n4e';
nodes = n4e(:,curElem);
coords = c4n(nodes,:);
coords = coords';
coords = reshape(coords,[2 3 length(curElem)]);
coords = permute(coords,[2 1 3]);
area = area4e(curElem)';

nrPoints = nrPts/length(curElem);
area = repmat(area,nrPoints,1);

P1X = repmat( squeeze(coords(1,1,:))',nrPoints,1 );
P2X = repmat( squeeze(coords(2,1,:))',nrPoints,1 );
P3X = repmat( squeeze(coords(3,1,:))',nrPoints,1 );

P1Y = repmat( squeeze(coords(1,2,:))',nrPoints,1 );
P2Y = repmat( squeeze(coords(2,2,:))',nrPoints,1 );
P3Y = repmat( squeeze(coords(3,2,:))',nrPoints,1 );

b1 = 1./2./area(:).*( (P2Y(:)-P3Y(:)).*x + (P3X(:)-P2X(:)).*y + P2X(:).*P3Y(:)-P3X(:).*P2Y(:) );
b2 = 1./2./area(:).*( (P3Y(:)-P1Y(:)).*x + (P1X(:)-P3X(:)).*y + P3X(:).*P1Y(:)-P1X(:).*P3Y(:) );
b3 = 1./2./area(:).*( (P1Y(:)-P2Y(:)).*x + (P2X(:)-P1X(:)).*y + P1X(:).*P2Y(:)-P2X(:).*P1Y(:) );

basisU = zeros(nrPts,2,6);
basisU(:,1,1) = b1;
basisU(:,1,2) = b2;
basisU(:,1,3) = b3;
basisU(:,2,4) = b1;
basisU(:,2,5) = b2;
basisU(:,2,6) = b3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basisP3 = getBasisP3(pts,curElem,lvl,p)
nrPts = size(pts,1);
basisU = p.statics.basisU;
basisP1 = basisU(pts,curElem,lvl,p);

b1 = basisP1(:,1,1)';
b2 = basisP1(:,1,2)';
b3 = basisP1(:,1,3)';

nrPoints = nrPts/length(curElem);
basisP3 = zeros(10,nrPoints,length(curElem));

basisP3(1,:,:) = b1;
basisP3(2,:,:) = b2;
basisP3(3,:,:) = b3;
basisP3(4,:,:) = b1.*b2;
basisP3(5,:,:) = b2.*b3;
basisP3(6,:,:) = b1.*b3;
basisP3(7,:,:) = b1.*b2.*(b1-b2);
basisP3(8,:,:) = b2.*b3.*(b2-b3);
basisP3(9,:,:) = b1.*b3.*(b3-b1);
basisP3(10,:,:) = b1.*b2.*b3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
