function p = ODP2getUhFuncs(p)
% author: David Guenther 

p.statics.u_h = @getU_h;
p.statics.grad_h = @getGrad_h;

%% supply u_h
function u_h = getU_h(points,curElem,lvl,p)

u = p.level(lvl).u4e;
basisU = p.statics.basisU;
basisU = basisU(points,curElem,lvl,p);

curU = u(curElem,:);
u_h = (curU * basisU')';

%% supply \grad u_h
function grad_h = getGrad_h(points,curElem,lvl,p)

stressBasis = p.statics.stressBasis;
u = p.level(lvl).u4e;

evalBasis = stressBasis(points,curElem,lvl,p);
u = ones(length(points(:,1)),1)*u(curElem,:);
u = reshape(u',[1 6 length(points(:,1))]);

grad_h = squeeze(matMul(u,evalBasis))';
