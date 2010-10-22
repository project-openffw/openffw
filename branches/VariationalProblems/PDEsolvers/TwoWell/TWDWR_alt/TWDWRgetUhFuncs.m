function p = TWDWRgetUhFuncs(p)
% author: Lena Noack

p.statics.u_h = @getU_h;
p.statics.Au_h = @getAu_h;

%p.statics.grad_h = @getGradU_h;
p.statics.gradAu_h = @getGradAu_h;

p.statics.I2u_h = @getI2u_h;
p.statics.gradI2u_h = @getGradI2u_h;
%% supply u_h
function u_h = getU_h(x,y,curElem,lvl,p)

u = p.level(lvl).u4e;
curU = u(curElem);
u_h = repmat(curU,length(x),1);

%% supply Au_{h}
function Auh = getAu_h(x,y,curElem,lvl,p)

Auh4n = p.level(lvl).Auh;
n4e = p.level(lvl).geom.n4e;
basisP1 = p.statics.basisP1;

evalBasisP1 = basisP1(x,y,curElem,lvl,p);

Auh = (Auh4n(n4e(curElem,:),1)'*evalBasisP1)';

%% supply \grad u
function grad_h = getGradU_h(x,y,curElem,lvl,p)
%No conj. Function defined!
%sigma_h = p.statics.sigma_h;

%CONV = p.params.CONV;
%if strcmp(CONV,'c')
%    NonLinearFuncDer1 = p.problem.nonLinearRegDerA;
%    NonLinearFuncDer2 = p.problem.nonLinearRegDerB;
%else
%    NonLinearFuncDer1 = p.problem.nonLinearExactDerA;
%    NonLinearFuncDer2 = p.problem.nonLinearExactDerB;
%end
%evalSigmah = sigma_h(x,y,curElem,lvl,p);
%absSigmah = (evalSigmah(:,1).^2 + evalSigmah(:,2).^2).^(1/2);
%evalNonLinear = NonLinearFuncDer1(absSigmah,curElem,lvl,p); %???
%grad_h = (evalNonLinear*[1,1]).*evalSigmah;

grad4e = p.level(lvl).grad4e;
curGrad = grad4e(curElem,:);

grad_h = zeros(length(x),2);
for j = 1:length(x)
    grad_h(j,:) = curGrad;
end

%% supply grad Au_{h}
function gradAuh = getGradAu_h(x,y,curElem,lvl,p)

Auh4n = p.level(lvl).Auh;
n4e = p.level(lvl).geom.n4e;

gradBasisP1 = p.statics.gradBasisP1;
evalGradBasisP1 = gradBasisP1(x,y,curElem,lvl,p);

Auh = Auh4n(n4e(curElem,:),1)';

gradAuh = squeeze(matMul(repmat(Auh,[1 1 length(x)]),evalGradBasisP1))';

%% supply I^2(\lambda2_h)
function I2Auh = getI2u_h(x,y,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];

Auh = p.statics.Au_h;
evalAuh = Auh(coords(:,1),coords(:,2),curElem,lvl,p);

basisP2 = p.statics.basisP2;
evalBasisP2 = basisP2(x,y,curElem,lvl,p);

I2Auh = (evalAuh' * evalBasisP2)';

%% supply grad I^2(\lambda2_h)
function gradI2uh = getGradI2u_h(x,y,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

Auh = p.statics.Au_h;
gradBasisP2 = p.statics.gradBasisP2;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];
      
evalAuh = Auh(coords(:,1),coords(:,2),curElem,lvl,p);

evalGradBasisP2 = gradBasisP2(x,y,curElem,lvl,p);

evalAuh = repmat(evalAuh',[1 1 length(x)]);

gradI2uh = squeeze(matMul(evalAuh,evalGradBasisP2));
