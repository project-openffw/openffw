function p = ODRTgetPhFuncs(p)
% author: David Guenther 
% Copyright 2007 David Guenther
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
%
%
%

p.statics.grad_h = @getGradU_h;
p.statics.sigma_h = @getSigma_h;
p.statics.Ap_h = @getAp_h;
p.statics.jacobianP_h = @getJacobianP_h;
p.statics.I2p_h = @getI2p_h;
p.statics.gradI2p_h = @getGradI2p_h;
p.statics.divI2p_h = @getDivI2p_h;

%% supply p_{\epsilon,h}
function sigma_h = getSigma_h(points,curElem,lvl,p)

stressBasis = p.statics.stressBasis;

grad4e = p.level(lvl).grad4e;
basisSigma = stressBasis(points,curElem,lvl,p);

sigma_h = zeros(length(points(:,1)),2);

for j = 1:length(points(:,1))
      sigma_h(j,:) = (grad4e(curElem,:)*basisSigma(:,:,j))';
end

%% supply \grad u_{\epsilon,h}
function grad_h = getGradU_h(points,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
conjNonLinearFuncDer = p.problem.conjNonLinearFuncDer;

evalSigmah = sigma_h(points,curElem,lvl,p);
absSigmah = (evalSigmah(:,1).^2 + evalSigmah(:,2).^2).^(1/2);

% \grad u = DW^*_\epsilon(p_h) = W^*'_\epsilon(|p_h|)/|p_h|*p_h
evalConjNonLinear = conjNonLinearFuncDer(absSigmah,curElem,lvl,p);

grad_h = (evalConjNonLinear*[1,1]).*evalSigmah;

%% supply Ap_{\epsilon,h}
function Aph = getAp_h(points,curElem,lvl,p)

basisP1 = p.statics.basisP1;
Aph4n = p.level(lvl).Aph;
n4e = p.level(lvl).geom.n4e;

evalBasisP1 = basisP1(points,curElem,lvl,p);

AphX = Aph4n(n4e(curElem,:),1)'*evalBasisP1';
AphY = Aph4n(n4e(curElem,:),2)'*evalBasisP1';

Aph = [AphX',AphY'];

%% supply jacobian of p_h
function jacobianP_h = getJacobianP_h(points,curElem,lvl,p)

grad4e = p.level(end).grad4e;
jacobianStressBasis = p.statics.jacobianStressBasis;
evalJacobian = jacobianStressBasis(points,curElem,lvl,p);

jacobianP_h = zeros(2,2,length(points(:,1)));

grad = grad4e(curElem,:);

for k = 1:length(points(:,1))
    jacobianP_h(:,:,k) = grad(1)*evalJacobian(:,:,1,k) + grad(2)*evalJacobian(:,:,2,k) +...
                         grad(3)*evalJacobian(:,:,3,k);
end

%% supply I^2(p_h)
function I2p_h = getI2p_h(points,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

Ap_h = p.statics.Ap_h;
basisP2 = p.statics.basisP2;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];
      
evalAph = Ap_h(coords(:,1),coords(:,2),curElem,lvl,p);

evalBasisP2 = basisP2(points,curElem,lvl,p);

I2p_h = (evalAph' * evalBasisP2)';

%% supply grad I^2(p_h)
function gradI2p_h = getGradI2p_h(points,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

Ap_h = p.statics.Ap_h;
gradBasisP2 = p.statics.gradBasisP2;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];
      
evalAph = Ap_h(coords(:,1),coords(:,2),curElem,lvl,p);

evalGradBasisP2 = gradBasisP2(points,curElem,lvl,p);

evalAph = repmat(evalAph',[1 1 length(points(:,1))]);

gradI2p_h = matMul(evalAph,evalGradBasisP2);

%% supply div I^2(p_h)
function divI2p_h = getDivI2p_h(points,curElem,lvl,p)

gradI2p_h = p.statics.gradI2p_h;

evalGradI2ph = gradI2p_h(points,curElem,lvl,p);

divI2p_h = squeeze(evalGradI2ph(1,1,:) + evalGradI2ph(2,2,:));
