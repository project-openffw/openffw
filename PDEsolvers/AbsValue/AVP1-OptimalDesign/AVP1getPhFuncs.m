function p = P1getPhFuncs(p)
%author: David Guenther

p.statics.sigma_h = @getSigma_h;
p.statics.Ap_h = @getAp_h;
p.statics.AvP_h = @getAvP_h;

%% supply p_h = DW(\grad u_h)
function sigma_h = getSigma_h(x,y,curElem,lvl,p)

nonLinearExactDer = p.problem.nonLinearExactDer;

grad4e = p.level(lvl).grad4e;
curGrad = grad4e(curElem,:);

sigma_h = nonLinearExactDer(norm(curGrad),curElem,lvl,p)*curGrad;
sigma_h = ones(length(x),1) * sigma_h;

%% supply average p_h: A(p_h), p_h=DW(grad Uh)
function Aph = getAp_h(x,y,curElem,lvl,p)

Aph4n = p.level(lvl).Aph;
n4e = p.level(lvl).geom.n4e;
basisU = p.statics.basisU;
basisP1 = basisU(x,y,curElem,lvl,p)';

AphX = Aph4n(n4e(curElem,:),1)'*basisP1;
AphY = Aph4n(n4e(curElem,:),2)'*basisP1;

Aph = [AphX',AphY'];

%% supply average p_h: A(p_h), p_h=grad Uh
function Aph = getAvP_h(x,y,curElem,lvl,p)

Aph4n = p.level(lvl).AvPh;
n4e = p.level(lvl).geom.n4e;
basisU = p.statics.basisU;
basisP1 = basisU(x,y,curElem,lvl,p)';

AphX = Aph4n(n4e(curElem,:),1)'*basisP1;
AphY = Aph4n(n4e(curElem,:),2)'*basisP1;

Aph = [AphX',AphY'];
