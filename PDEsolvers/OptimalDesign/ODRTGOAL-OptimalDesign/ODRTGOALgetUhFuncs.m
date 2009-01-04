function p = ODRTGOALgetUhFuncs(p)
% author: David Guenther

p.statics.u_h = @getU_h;
p.statics.Au_h = @getAu_h;

p.statics.grad_h = @getGradU_h;
p.statics.gradAu_h = @getGradAu_h;

p.statics.I2u_h = @getI2u_h;
p.statics.gradI2u_h = @getGradI2u_h;
%% supply u_h
function u_h = getU_h(points,curElem,lvl,p)

u = p.level(lvl).u4e;
curU = u(curElem);
u_h = repmat(curU,length(points(:,1)),1);

%% supply Au_{h}
function Auh = getAu_h(points,curElem,lvl,p)

Auh4n = p.level(lvl).Auh;
n4e = p.level(lvl).geom.n4e;
basisP1 = p.statics.basisP1;

evalBasisP1 = basisP1(points,curElem,lvl,p);

Auh = (Auh4n(n4e(curElem,:),1)'*evalBasisP1)';

%% supply \grad u_{\epsilon,h}
function grad_h = getGradU_h(points,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
conjNonLinearFuncDer = p.problem.conjNonLinearFuncDer;

evalSigmah = sigma_h(points,curElem,lvl,p);
absSigmah = (evalSigmah(:,1).^2 + evalSigmah(:,2).^2).^(1/2);

evalConjNonLinear = conjNonLinearFuncDer(absSigmah,curElem,lvl,p);

grad_h = (evalConjNonLinear*[1,1]).*evalSigmah;

%% supply grad Au_{h}
function gradAuh = getGradAu_h(points,curElem,lvl,p)

Auh4n = p.level(lvl).Auh;
n4e = p.level(lvl).geom.n4e;

gradBasisP1 = p.statics.gradBasisP1;
evalGradBasisP1 = gradBasisP1(points,curElem,lvl,p);

Auh = Auh4n(n4e(curElem,:),1)';

gradAuh = squeeze(matMul(repmat(Auh,[1 1 length(points(:,1))]),evalGradBasisP1))';

%% supply I^2(\lambda2_h)
function I2Auh = getI2u_h(points,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];

Auh = p.statics.Au_h;
evalAuh = Auh(coords(:,1),coords(:,2),curElem,lvl,p);

basisP2 = p.statics.basisP2;
evalBasisP2 = basisP2(points,curElem,lvl,p);

I2Auh = (evalAuh' * evalBasisP2)';

%% supply grad I^2(\lambda2_h)
function gradI2uh = getGradI2u_h(points,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

Auh = p.statics.Au_h;
gradBasisP2 = p.statics.gradBasisP2;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];
      
evalAuh = Auh(coords(:,1),coords(:,2),curElem,lvl,p);

evalGradBasisP2 = gradBasisP2(points,curElem,lvl,p);

evalAuh = repmat(evalAuh',[1 1 length(points(:,1))]);

gradI2uh = squeeze(matMul(evalAuh,evalGradBasisP2));
