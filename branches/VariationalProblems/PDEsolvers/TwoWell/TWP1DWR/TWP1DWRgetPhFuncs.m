function p = P1getPhFuncs(p)
%author: David Guenther, Lena Noack

p.statics.sigma_h = @getSigma_h;
p.statics.Ap_h = @getAp_h;

%% supply p_h = DW(\grad u_h)
function sigma_h = getSigma_h(x,y,curElem,lvl,p)
F1 = p.problem.F1;
F2 = p.problem.F2;

CONV = p.params.CONV;
if strcmp(CONV,'c')
nonLinearRegDer1 = p.problem.nonLinearRegDerA; %W'**1(F)
nonLinearRegDer2 = p.problem.nonLinearRegDerB; %W'**2(F)
else
nonLinearExactDer1 = p.problem.nonLinearExactDerA; %W'1(F)
nonLinearExactDer2 = p.problem.nonLinearExactDerB; %W'2(F)
end

grad4e = p.level(lvl).grad4e;
curGrad = grad4e(curElem,:);

if strcmp(CONV,'c')
%% sigma_h = DW**(F) = W'**1(F)F + W'**2(F)(F2,F)F2
sigma_h = nonLinearRegDer1(curGrad(1),curGrad(2),curElem,lvl,p)*curGrad +...
    nonLinearRegDer2(curGrad(1),curGrad(2),curElem,lvl,p)*(curGrad(1)*F2(1)+curGrad(2)*F2(2))*F2';
else
%% sigma_h = DW(F) = W'1(F)(F-F1) + W'2(F)(F-F2)
sigma_h = nonLinearExactDer1(curGrad(1),curGrad(2),curElem,lvl,p)*(curGrad-F1') +...
    nonLinearExactDer2(curGrad(1),curGrad(2),curElem,lvl,p)*(curGrad-F2');
end

sigma_h = ones(length(x),1) * sigma_h;

%% supply average p_h: A(p_h)
function Aph = getAp_h(x,y,curElem,lvl,p)

Aph4n = p.level(lvl).Aph;
n4e = p.level(lvl).geom.n4e;
basisU = p.statics.basisU;
basisP1 = basisU(x,y,curElem,lvl,p)';

AphX = Aph4n(n4e(curElem,:),1)'*basisP1;
AphY = Aph4n(n4e(curElem,:),2)'*basisP1;

Aph = [AphX',AphY'];

