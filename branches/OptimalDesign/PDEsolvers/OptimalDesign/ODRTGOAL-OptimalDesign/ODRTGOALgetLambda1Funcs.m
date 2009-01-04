function p = ODRTGOALgetLambda1Funcs(p)
% author: David Guenther

p.statics.lambda1 = @getLambda1;
p.statics.divLambda1_h = @getDivLambda1_h;
p.statics.divAlambda1_h = @getDivAlambda1;
p.statics.Alambda1_h = @getAlambda1_h;
p.statics.I2lambda1_h = @getI2lambda1_h;
p.statics.gradI2lambda1_h = @getGradI2lambda1_h;
p.statics.divI2lambda1_h = @getDivI2lambda1_h;

%% supply \lambda1
function lambda_h = getLambda1(points,curElem,lvl,p)

lambda14e = p.level(lvl).lambda14e;
basisSigma = p.statics.stressBasis;
evalBasisSigma = basisSigma(points,curElem,lvl,p);

lambda_h = zeros(length(points(:,1)),2);

for j = 1:length(points(:,1))
      lambda_h(j,:) = (lambda14e(curElem,:)*evalBasisSigma(:,:,j))';
end

%% supply div(\lambda1_h)
function divLambda1 = getDivLambda1_h(points,curElem,lvl,p)

divStressBasis = p.statics.divStressBasis;
lambda14e = p.level(lvl).grad4e;

evalBasis = divStressBasis(points,curElem,lvl,p);

divLambda1 = (lambda14e(curElem,:)*evalBasis)';

%% supply Alambda1
function Alambda1 = getAlambda1_h(points,curElem,lvl,p)

Alambda14n = p.level(lvl).Alambda1;
n4e = p.level(lvl).geom.n4e;

basisP1 = p.statics.basisP1;
evalBasisP1 = basisP1(points,curElem,lvl,p);

Alambda1X = Alambda14n(n4e(curElem,:),1)'*evalBasisP1;
Alambda1Y = Alambda14n(n4e(curElem,:),2)'*evalBasisP1;

Alambda1 = [Alambda1X',Alambda1Y'];

%% supply divAlambda_1
function divAlambda1 = getDivAlambda1(points,curElem,lvl,p)

Alambda14n = p.level(lvl).Alambda1;
n4e = p.level(lvl).geom.n4e;

gradBasisP1 = p.level(lvl).enum.grad4e;

Alambda1X = Alambda14n(n4e(curElem,:),1)'*gradBasisP1(:,:,curElem);
Alambda1Y = Alambda14n(n4e(curElem,:),2)'*gradBasisP1(:,:,curElem);

divAlambda1 = Alambda1X(1) + Alambda1Y(2);
divAlambda1 = divAlambda1*ones(length(points(:,1)),1);

%% supply I^2(\lambda1_h)
function I2lambda1 = getI2lambda1_h(points,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

Alambda1 = p.statics.Alambda1;
basisP2 = p.statics.basisP2;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];
      
evalAlambda1 = Alambda1(coords,curElem,lvl,p);

evalBasisP2 = basisP2(points,curElem,lvl,p);

I2lambda1 = (evalAlambda1' * evalBasisP2)';

%% supply grad I^2(\lambda1_h)
function gradI2lambda1 = getGradI2lambda1_h(points,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

Alambda1 = p.statics.Alambda1;
gradBasisP2 = p.statics.gradBasisP2;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];
      
evalAlambda1 = Alambda1(coords,curElem,lvl,p);

evalGradBasisP2 = gradBasisP2(points,curElem,lvl,p);

evalAlambda1 = repmat(evalAlambda1',[1 1 length(points(:,1))]);

gradI2lambda1 = matMul(evalAlambda1,evalGradBasisP2);

%% supply div I^2(\lambda1_h)
function divI2lambda1_h = getDivI2lambda1_h(points,curElem,lvl,p)

gradI2lambda1_h = p.statics.gradI2lambda1_h;

evalGradI2lambda1h = gradI2lambda1_h(points,curElem,lvl,p);

divI2lambda1_h = squeeze(evalGradI2lambda1h(1,1,:) + evalGradI2lambda1h(2,2,:));