function p = elliptic_init(p)
% initialisations for elliptic pde's

% Copyright 2007 Jan Reininghaus, David Guenther, 
%                Andreas Byfut, Joscha Gedicke
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
p.statics.H1semiError = @H1semiErrorElliptic;
p.statics.L2error = @L2errorElliptic;
p.statics.sigma_h = @getSigma_h;
p.statics.u_h = @getU_h;
p.statics.gradU_h = @getGradU_h;
p.statics.d2U_h = @getD2U_h;

p.statics.H1semiErrorVectorised = @H1semiErrorEllipticVectorised;
p.statics.L2errorVectorised = @L2errorEllipticVectorised;
p.statics.sigma_hVectorised = @getSigma_hVectorised;
p.statics.u_hVectorised = @getU_hVectorised;
p.statics.gradU_hVectorised = @getGradU_hVectorised;
p.statics.d2U_hVectorised = @getD2U_hVectorised;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = L2errorElliptic(pts,curElem,lvl,p)
u_h = p.statics.u_h;
u_exact = p.problem.u_exact;
exactU =  u_exact(pts,p)';
uh = u_h(pts,curElem,lvl,p);
val(1,1,:) = (uh(:) - exactU(:)).*(uh(:) - exactU(:));


%%
function val = H1semiErrorElliptic(pts,curElem,lvl,p)
gradU_exact = p.problem.gradU_exact;
kappa = p.problem.kappa;
curKappa = kappa(pts,p);
curGradU_exact =  gradU_exact(pts,p);
curGradU_exact = reshape (curGradU_exact',2,1,[]);
exactSigma = matMul(curKappa,curGradU_exact);
exactSigma = reshape(exactSigma,2,[])';
sigma_h = p.statics.sigma_h(pts,curElem,curKappa,lvl,p);
val(1,1,:) = sum( (sigma_h - exactSigma).*(sigma_h - exactSigma),2 );

%%
function val = getSigma_h(pts,curElem,curKappa,lvl,p)
gradU_h = p.statics.gradU_h;
curGrad = permute(gradU_h(pts,curElem,lvl,p),[2 1 3]);
sigma_h = matMul(curKappa,curGrad);
sigma_h  = reshape(sigma_h,2,[])';
val = sigma_h;

%%
function val = getU_h(pts,curElem,lvl,p)
curU = p.level(lvl).u4e(curElem,:);
basis = p.statics.basis(pts,curElem,lvl,p);
val(1,1,:) = curU * basis';

%%
function val = getGradU_h(pts,curElem,lvl,p)
nrPts = size(pts,1);
curU = p.level(lvl).u4e(curElem,:);
curGrad = p.statics.gradBasis(pts,curElem,lvl,p);
curU = reshape(curU(:)*ones(1,nrPts),[size(curU,1),size(curU,2),nrPts]);
val = zeros([1 2 nrPts]);
val(1,:,:) = matMul(curU,curGrad);

%%
function val = getD2U_h(pts,curElem,lvl,p)
nrPts = size(pts,1);
curU = p.level(lvl).u4e(curElem,:);
curD2 = p.statics.d2Basis(pts,curElem,lvl,p);
curU = reshape(curU'*ones(1,nrPts),[1 size(curU,2) nrPts]);
curD2_1(:,1,:) = curD2(1,:,:);
curD2_2(:,1,:) = curD2(2,:,:);
curD2_3(:,1,:) = curD2(3,:,:);
val = zeros(2,2,nrPts);
val(1,1,:) = matMul(curU,curD2_1);
val(2,2,:) = matMul(curU,curD2_2);
val(1,2,:) = matMul(curU,curD2_3);
val(2,1,:) = val(1,2,:);

%% Vectorised %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function val = L2errorEllipticVectorised(pts,pts_ref,parts,lvl,p)
u_h = p.statics.u_hVectorised;
u_exact = p.problem.u_exact;
exactU =  u_exact(pts,p);
uh =  u_h(pts,pts_ref,parts,lvl,p);
val(1,1,:) = (uh(:) - exactU(:)).*(uh(:) - exactU(:));


%%
function val = H1semiErrorEllipticVectorised(pts,pts_ref,parts,lvl,p)
gradU_h = p.statics.gradU_hVectorised;
gradU_exact = p.problem.gradU_exact;
kappa = p.problem.kappa;
curKappa = kappa(pts,p);
curGradU_exact =  gradU_exact(pts,p);
curGradU_exact = reshape (curGradU_exact',2,1,[]);
exactSigma = matMul(curKappa,curGradU_exact);
exactSigma = reshape(exactSigma,2,[])';
curGrad = permute(gradU_h(pts,pts_ref,parts,lvl,p),[2 1 3]);
sigma_h = matMul(curKappa,curGrad);
sigma_h = reshape(sigma_h,2,[])';
val(1,1,:) = sum( (sigma_h - exactSigma).*(sigma_h - exactSigma),2 );

%%
function val = getSigma_hVectorised(pts,pts_ref,parts,lvl,p)
kappa = p.problem.kappa;
gradU_h = p.statics.gradU_hVectorised;
curKappa = kappa(pts,p);
curGrad = permute(gradU_h(pts,pts_ref,parts,lvl,p),[2 1 3]);
sigma_h = matMul(curKappa,curGrad);
sigma_h = reshape(sigma_h,2,[]);
val(1,:,:) = sigma_h;

%%
function val = getU_hVectorised(pts,pts_ref,parts,lvl,p)
u = p.level(lvl).u4e(parts,:);
basis = p.statics.basisVectorised(pts,pts_ref,parts,lvl,p);
val(1,1,:) = sum(u .* basis,2);

%%
function val = getGradU_hVectorised(pts,pts_ref,parts,lvl,p)
u = p.level(lvl).u4e(parts,:);
gradBasis = p.statics.gradBasisVectorised(pts,pts_ref,parts,lvl,p);
u = reshape(u',[1 size(u,2) size(u,1) ]);
val(1,:,:) = sum(matMul(u,gradBasis),1);

%%
function val = getD2U_hVectorised(pts,pts_ref,parts,lvl,p)
nrPts = size(pts,1);
u =  p.level(lvl).u4e(parts,:);
d2Basis = p.statics.d2BasisVectorised(pts,pts_ref,parts,lvl,p);
u = reshape(u',[1 size(u,2) size(u,1) ]);
d2Basis1 = permute(d2Basis(1,:,:),[2 1 3]);
d2Basis2 = permute(d2Basis(2,:,:),[2 1 3]);
d2Basis3 = permute(d2Basis(3,:,:),[2 1 3]);
val = zeros(2,2,nrPts);
val(1,1,:) = matMul(u,d2Basis1);
val(2,2,:) = matMul(u,d2Basis2);
val(1,2,:) = matMul(u,d2Basis3);
val(2,1,:) = val(1,2,:);