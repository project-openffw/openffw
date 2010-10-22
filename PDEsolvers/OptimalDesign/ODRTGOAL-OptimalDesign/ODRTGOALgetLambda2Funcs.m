function p = ODRTGOALgetLambda2Funcs(p)
% author: David Guenther

p.statics.lambda2 = @getLambda2;
p.statics.Alambda2_h = @getAlambda2;
p.statics.I2lambda2_h = @getI2lambda2_h;
p.statics.gradI2lambda2_h = @getGradI2lambda2_h;

%% supply \lambda2
function lambda_h = getLambda2(x,y,curElem,lvl,p)

u = p.level(lvl).lambda24e;
curU = u(curElem);
lambda_h = repmat(curU,length(x),1);

%% supply A\lambda2
function Alambda2 = getAlambda2(x,y,curElem,lvl,p)

Alambda24n = p.level(lvl).Alambda2;
n4e = p.level(lvl).geom.n4e;
basisP1 = p.statics.basisP1;

evalBasisP1 = basisP1(x,y,curElem,lvl,p);

Alambda2 = (Alambda24n(n4e(curElem,:))'*evalBasisP1)';

%% supply I^2(\lambda2_h)
function I2lambda2 = getI2lambda2_h(x,y,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];

Alambda2 = p.statics.Alambda2;
evalAlambda2 = Alambda2(coords(:,1),coords(:,2),curElem,lvl,p);

basisP2 = p.statics.basisP2;
evalBasisP2 = basisP2(x,y,curElem,lvl,p);

I2lambda2 = (evalAlambda2' * evalBasisP2)';

%% supply grad I^2(\lambda2_h)
function gradI2lambda2 = getGradI2lambda2_h(x,y,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

Alambda2 = p.statics.Alambda2;
gradBasisP2 = p.statics.gradBasisP2;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];
      
evalAlambda2 = Alambda2(coords(:,1),coords(:,2),curElem,lvl,p);

evalGradBasisP2 = gradBasisP2(x,y,curElem,lvl,p);

evalAlambda2 = repmat(evalAlambda2',[1 1 length(x)]);

gradI2lambda2 = squeeze(matMul(evalAlambda2,evalGradBasisP2));
