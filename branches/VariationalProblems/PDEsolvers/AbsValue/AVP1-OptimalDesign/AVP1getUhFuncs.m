function p = P1getUhFuncs(p)
%author: David Guenther

p.statics.u_h = @getU_h;
p.statics.grad_h = @getGrad_h;
p.statics.AGrad_h = @getAGrad_h;

%% supply u_h
function u_h = getU_h(x,y,curElem,lvl,p)

u = p.level(lvl).u4e;
basisU = p.statics.basisU;
basisU = basisU(x,y,curElem,lvl,p)';

curU = u(curElem,:);
u_h = (curU * basisU)';

%% supply \grad u_h
function grad_h = getGrad_h(x,y,curElem,lvl,p)

grad4e = p.level(lvl).grad4e;
curGrad = grad4e(curElem,:);

grad_h = zeros(length(x),2);
for j = 1:length(x)
    grad_h(j,:) = curGrad;
end

%% supply average p_h: A(p_h)
function Aph = getAGrad_h(x,y,curElem,lvl,p)

Aph4n = p.level(lvl).AGradh;
n4e = p.level(lvl).geom.n4e;
basisU = p.statics.basisU;
basisP1 = basisU(x,y,curElem,lvl,p)';

AphX = Aph4n(n4e(curElem,:),1)'*basisP1;
AphY = Aph4n(n4e(curElem,:),2)'*basisP1;

Aph = [AphX',AphY'];
