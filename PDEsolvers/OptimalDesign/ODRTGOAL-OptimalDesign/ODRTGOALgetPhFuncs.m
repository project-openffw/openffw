function p = ODRTGOALgetPhFuncs(p)
% author: David Guenther

p.statics.sigma_h = @getSigma_h;
p.statics.divP_h = @getDivP_h;
p.statics.Ap_h = @getAp_h;
p.statics.divAp_h = @getDivAp_h;
p.statics.I2p_h = @getI2p_h;
p.statics.gradI2p_h = @getGradI2p_h;
p.statics.divI2p_h = @getDivI2p_h;

%% supply p_{\epsilon,h}
function sigma_h = getSigma_h(points,curElem,lvl,p)

grad4e = p.level(lvl).grad4e;

basisSigma = p.statics.stressBasis;
evalBasisSigma = basisSigma(points,curElem,lvl,p);

sigma_h = zeros(length(points(:,1)),2);

for j = 1:length(points(:,1))
      sigma_h(j,:) = (grad4e(curElem,:)*evalBasisSigma(:,:,j))';
end

%% supply div(p_h)
function divPh = getDivP_h(points,curElem,lvl,p)

divStressBasis = p.statics.divStressBasis;
grad4e = p.level(lvl).grad4e;

evalBasis = divStressBasis(points,curElem,lvl,p);

divPh = (grad4e(curElem,:)*evalBasis)';

%% supply Ap_{\epsilon,h}
function Aph = getAp_h(points,curElem,lvl,p)

Aph4n = p.level(lvl).Aph;
n4e = p.level(lvl).geom.n4e;

basisP1 = p.statics.basisP1;
evalBasisP1 = basisP1(points,curElem,lvl,p);

AphX = Aph4n(n4e(curElem,:),1)'*evalBasisP1;
AphY = Aph4n(n4e(curElem,:),2)'*evalBasisP1;

Aph = [AphX',AphY'];

%% supply divAp_{\epsilon,h}
function divAph = getDivAp_h(points,curElem,lvl,p)

Aph4n = p.level(lvl).Aph;
n4e = p.level(lvl).geom.n4e;

gradBasisP1 = p.level(lvl).enum.grad4e;

AphX = Aph4n(n4e(curElem,:),1)'*gradBasisP1(:,:,curElem);
AphY = Aph4n(n4e(curElem,:),2)'*gradBasisP1(:,:,curElem);

divAph = AphX(1) + AphY(2);
divAph = divAph *ones(length(points(:,1)),1);

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
      
evalAph = Ap_h(coords,curElem,lvl,p);

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
      
evalAph = Ap_h(coords,curElem,lvl,p);

evalGradBasisP2 = gradBasisP2(points,curElem,lvl,p);

evalAph = repmat(evalAph',[1 1 length(points(:,1))]);

gradI2p_h = matMul(evalAph,evalGradBasisP2);

%% supply div I^2(p_h)
function divI2p_h = getDivI2p_h(points,curElem,lvl,p)

gradI2p_h = p.statics.gradI2p_h;

evalGradI2ph = gradI2p_h(points,curElem,lvl,p);

divI2p_h = squeeze(evalGradI2ph(1,1,:) + evalGradI2ph(2,2,:));