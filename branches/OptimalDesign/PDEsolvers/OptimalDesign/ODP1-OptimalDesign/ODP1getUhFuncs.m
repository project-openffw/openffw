function p = ODP1getUhFuncs(p)
%author: David Guenther

p.statics.u_h = @getU_h;
p.statics.grad_h = @getGrad_h;

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
