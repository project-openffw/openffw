function p = ODP1getUhFuncs(p)
%author: David Guenther

p.statics.u_h = @getU_h;
p.statics.grad_h = @getGrad_h;

%% supply u_h
function u_h = getU_h(points,curElem,lvl,p)

u = p.level(lvl).u4e;
basisU = p.statics.basisU;

basisU = basisU(points,curElem,lvl,p)';

curU = u(curElem,:);
u_h = (curU * basisU)';

%% supply \grad u_h
function grad_h = getGrad_h(points,curElem,lvl,p)

grad4e = p.level(lvl).grad4e;
curGrad = grad4e(curElem,:);

grad_h = zeros(length(points(:,1)),2);
for j = 1:length(points(:,1))
    grad_h(j,:) = curGrad;
end
