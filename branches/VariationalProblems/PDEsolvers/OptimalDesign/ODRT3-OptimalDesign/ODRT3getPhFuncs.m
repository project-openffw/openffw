function p = ODRTgetPhFuncs(p)
% author: David Guenther 

p.statics.grad_h = @getGradU_h;
p.statics.DWRgrad_h = @getDWRGradU_h;
p.statics.sigma_h = @getSigma_h;
p.statics.DWRsigma_h = @getDWRSigma_h;
p.statics.Ap_h = @getAp_h;
p.statics.DWRAp_h = @getDWRAp_h;
p.statics.jacobianP_h = @getJacobianP_h;
p.statics.I2p_h = @getI2p_h;
p.statics.gradI2p_h = @getGradI2p_h;
p.statics.divI2p_h = @getDivI2p_h;

%% supply p_{\epsilon,h}
function sigma_h = getSigma_h(x,y,curElem,lvl,p)

stressBasis = p.statics.stressBasis;

grad4e = p.level(lvl).grad4e;
basisSigma = stressBasis(x,y,curElem,lvl,p);

sigma_h = zeros(length(x),2);

for j = 1:length(x)
      sigma_h(j,:) = (grad4e(curElem,:)*basisSigma(:,:,j))';
end

%% supply dual p_{\epsilon,h}
function sigma_h = getDWRSigma_h(x,y,curElem,lvl,p)

stressBasis = p.statics.stressBasis;

grad4e = p.level(lvl).DWRgrad4e;
basisSigma = stressBasis(x,y,curElem,lvl,p);

sigma_h = zeros(length(x),2);

for j = 1:length(x)
      sigma_h(j,:) = (grad4e(curElem,:)*basisSigma(:,:,j))';
end

%% supply \grad u_{\epsilon,h}
function grad_h = getGradU_h(x,y,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
conjNonLinearFuncDer = p.problem.conjNonLinearFuncDer;

evalSigmah = sigma_h(x,y,curElem,lvl,p);
absSigmah = (evalSigmah(:,1).^2 + evalSigmah(:,2).^2).^(1/2);

% \grad u = DW^*_\epsilon(p_h) = W^*'_\epsilon(|p_h|)/|p_h|*p_h
%evalConjNonLinear = conjNonLinearFuncDer(evalSigmah(:,1),evalSigmah(:,2),curElem,lvl,p);
evalConjNonLinear = conjNonLinearFuncDer(absSigmah,curElem,lvl,p);

grad_h = (evalConjNonLinear*[1,1]).*evalSigmah;

%% supply dual \grad u_{\epsilon,h}
function DWRgrad_h = getDWRGradU_h(x,y,curElem,lvl,p)
% for dual solutions z,tau the relation Dz=DW^*(tau) cannot be garanteared

%sigma_h = p.statics.DWRsigma_h;
%conjNonLinearFuncDer = p.problem.conjNonLinearFuncDer;
%evalSigmah = sigma_h(x,y,curElem,lvl,p);
%absSigmah = (evalSigmah(:,1).^2 + evalSigmah(:,2).^2).^(1/2);
%evalConjNonLinear = conjNonLinearFuncDer(absSigmah,curElem,lvl,p);
%DWRgrad_h = (evalConjNonLinear*[1,1]).*evalSigmah;

basisU = p.statics.stressBasis;
basisU = basisU(x,y,curElem,lvl,p);

grad4e = p.level(lvl).DWRgradU4e;
curGrad = grad4e(curElem,:);

DWRgrad_h = zeros(length(x),2);
for j = 1:length(x)
    DWRgrad_h(j,:) = curGrad * basisU(:,:,j);
end

%% supply Ap_{\epsilon,h}
function Aph = getAp_h(x,y,curElem,lvl,p)

basisP1 = p.statics.basisP1;
Aph4n = p.level(lvl).Aph;
n4e = p.level(lvl).geom.n4e;

evalBasisP1 = basisP1(x,y,curElem,lvl,p);

AphX = Aph4n(n4e(curElem,:),1)'*evalBasisP1';
AphY = Aph4n(n4e(curElem,:),2)'*evalBasisP1';

Aph = [AphX',AphY'];

%% supply Ap_{\epsilon,h}
function Aph = getDWRAp_h(x,y,curElem,lvl,p)

basisP1 = p.statics.basisP1;
Aph4n = p.level(lvl).DWRAph;
n4e = p.level(lvl).geom.n4e;

evalBasisP1 = basisP1(x,y,curElem,lvl,p);

AphX = Aph4n(n4e(curElem,:),1)'*evalBasisP1';
AphY = Aph4n(n4e(curElem,:),2)'*evalBasisP1';

Aph = [AphX',AphY'];

%% supply jacobian of p_h
function jacobianP_h = getJacobianP_h(x,y,curElem,lvl,p)

grad4e = p.level(end).grad4e;
jacobianStressBasis = p.statics.jacobianStressBasis;
evalJacobian = jacobianStressBasis(x,y,curElem,lvl,p);

jacobianP_h = zeros(2,2,length(x));

grad = grad4e(curElem,:);

for k = 1:length(x)
    jacobianP_h(:,:,k) = grad(1)*evalJacobian(:,:,1,k) + grad(2)*evalJacobian(:,:,2,k) +...
                         grad(3)*evalJacobian(:,:,3,k);
end

%% supply I^2(p_h)
function I2p_h = getI2p_h(x,y,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

Ap_h = p.statics.Ap_h;
basisP2 = p.statics.basisP2;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];
      
evalAph = Ap_h(coords(:,1),coords(:,2),curElem,lvl,p);

evalBasisP2 = basisP2(x,y,curElem,lvl,p);

I2p_h = (evalAph' * evalBasisP2)';

%% supply grad I^2(p_h)
function gradI2p_h = getGradI2p_h(x,y,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

Ap_h = p.statics.Ap_h;
gradBasisP2 = p.statics.gradBasisP2;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];
      
evalAph = Ap_h(coords(:,1),coords(:,2),curElem,lvl,p);

evalGradBasisP2 = gradBasisP2(x,y,curElem,lvl,p);

evalAph = repmat(evalAph',[1 1 length(x)]);

gradI2p_h = matMul(evalAph,evalGradBasisP2);

%% supply div I^2(p_h)
function divI2p_h = getDivI2p_h(x,y,curElem,lvl,p)

gradI2p_h = p.statics.gradI2p_h;

evalGradI2ph = gradI2p_h(x,y,curElem,lvl,p);

divI2p_h = squeeze(evalGradI2ph(1,1,:) + evalGradI2ph(2,2,:));
