function p = ODRTGOALestimate(p)
% author: David Guenther
%% Input
nrElems = p.level(end).nrElems;
lvl = size(p.level,2);
n4e = p.level(end).geom.n4e;

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% supply error estimator
rho1Dual = zeros(nrElems,1);
rho2Dual = zeros(nrElems,1);
rho1 = zeros(nrElems,1);
rho2 = zeros(nrElems,1);
dummy1 = zeros(nrElems,1);
dummy2 = zeros(nrElems,1);

for curElem = 1:nrElems
    %compute DW(p)(Ap-p)
    DW_1 = integrate(n4e(curElem,:),lvl,degree,@integrandDW1,p);
    %compute D2W(p)(lambda1,Ap-p)
    D2W = integrate(n4e(curElem,:),lvl,degree,@integrandD2W,p);
    %compute lambda2*div(Ap-p)
    term1 = integrate(n4e(curElem,:),lvl,degree,@integrand6,p);
    
    dummy1(curElem) = term1;
    %set rho1*(p,u)(Ap-p)
    rho1Dual(curElem) = DW_1 - D2W - term1;
    
    %set rho2*(p,u)(Au-u)
    rho2Dual(curElem) = -integrate(n4e(curElem,:),lvl,degree,@integrand2,p); 
    %compute DW(p,u)(Alambda1-lambda1)
    DW_2 = integrate(n4e(curElem,:),lvl,degree,@integrandDW2,p);
    %compute u*div(Alambda1-lambda1)
    term2 = integrate(n4e(curElem,:),lvl,degree,@integrand4,p);
    
    dummy2(curElem) = term2;    
    %set rho1(p,u)(Alambda1-lambda1)
    rho1(curElem) = -DW_2 - term2;
    
     %compute f*(Alambda2-lambda2)
    term3 = integrate(n4e(curElem,:),lvl,degree,@integrand7,p);    
    %set rho2(p,u)(Alambda2-lambda2)
    rho2(curElem) = term3 - integrate(n4e(curElem,:),lvl,degree,@integrand5,p); 
end

residual = (1/2*(rho1Dual + rho2Dual) + 1/2*(rho1 + rho2)).^2;
%linError = (1/2*((rho1Dual + rho2Dual) - (rho1 + rho2))).^2;
   %-> eta = 0.5(rho*_1 + rho*_2 + rho_1 + rho_2) + Delta rho
linError = 0;%-> eta = 0.5(rho*_1 + rho*_2 + rho_1 + rho_2)

etaT = sqrt(residual);
etaOsc = sqrt(linError); 

%% Output
p.level(end).rho1Dual = rho1Dual;
p.level(end).rho2Dual = rho2Dual;
p.level(end).rho1 = rho1;
p.level(end).rho2 = rho2;
p.level(end).etaT = etaT;
p.level(end).etaOsc = etaOsc;
p.level(end).residual = sqrt(residual);
p.level(end).linError = sqrt(linError);
p.level(end).estimatedError = sqrt(sum(residual) + sum(linError));

%% supply integrand: DW_1
function val = integrandDW1(x,y,curElem,lvl,p)

grad_h = p.statics.grad_h;
evalIntegrand1 = integrand1(x,y,curElem,lvl,p);
evalGrad = grad_h(x,y,curElem,lvl,p);

val = sum(evalGrad.*evalIntegrand1,2);
val = reshape(val,[1 1 length(x)]);

%% supply integrand: DW_2
function val = integrandDW2(x,y,curElem,lvl,p)

grad_h = p.statics.grad_h;
evalIntegrand3 = integrand3(x,y,curElem,lvl,p);
evalGrad = grad_h(x,y,curElem,lvl,p);

val = sum(evalGrad.*evalIntegrand3,2);
val = reshape(val,[1 1 length(x)]);

%% supply integrand
function val = integrandD2W(x,y,curElem,lvl,p)
% W''(|X|)/|X|^2*(X*Y)(X*Z) + W'(|X|)/|X|^3*( Y*Z*|X|^2 - (X*Y)(X*Z) )

DW = p.problem.conjNonLinearFuncDer;
D2W = p.problem.conjNonLinearFuncSecDer;
sigma_h = p.statics.sigma_h;
lambda1 = p.statics.lambda1;


evalSigma = sigma_h(x,y,curElem,lvl,p);
evalLambda = lambda1(x,y,curElem,lvl,p);
evalIntegrand1 = integrand1(x,y,curElem,lvl,p);

absSigma = ( evalSigma(:,1).^2 + evalSigma(:,2).^2 ).^(1/2);
evalDW = DW(absSigma,curElem,lvl,p);
evalD2W = D2W(absSigma,curElem,lvl,p);

evalSigma = reshape(evalSigma',[1 2 length(x)]);
evalLambda = reshape(evalLambda',[1 2 length(x)]);
evalIntegrand1 = reshape(evalIntegrand1',[1 2 length(x)]);

YZ = matMul(evalLambda,permute(evalIntegrand1,[2 1 3]));
XY = matMul(evalSigma,permute(evalLambda,[2 1 3]));
XZ = matMul(evalSigma,permute(evalIntegrand1,[2 1 3]));
XYXZ = matMul(permute(XY,[2 1 3]),XZ);

if norm(absSigma) > 0
    term1 = matMul(reshape(evalD2W./absSigma.^2,[1 1 length(x)]),XYXZ);
    term2 = -matMul(reshape(evalDW./absSigma.^2,[1 1 length(x)]),XYXZ);
else
    term1 = 0;
    term2 = 0;
end

term3 = matMul(reshape(evalDW,[1 1 length(x)]),YZ);

val = term1 + term2 + term3;

%% supply first integrand: Ap-p
function val = integrand1(x,y,curElem,lvl,p)

p_h = p.statics.sigma_h;
Ap_h = p.statics.Ap_h;

evalPh = p_h(x,y,curElem,lvl,p);
evalAph = Ap_h(x,y,curElem,lvl,p);

val = evalAph - evalPh;

%% supply second integrand: div(lambda1)(Au-u)
function val = integrand2(x,y,curElem,lvl,p)

divLambda1_h = p.statics.divLambda1_h;
u_h = p.statics.u_h;
Au_h = p.statics.Au_h;

evalDivLambda1 = divLambda1_h(x,y,curElem,lvl,p);
evalU = u_h(x,y,curElem,lvl,p);
evalAu = Au_h(x,y,curElem,lvl,p);

val = evalDivLambda1.*(evalAu - evalU);
val = reshape(val,[1,1,length(x)]);

%% supply third integrand: A\lambda1-lambda1
function val = integrand3(x,y,curElem,lvl,p)

lambda1_h = p.statics.lambda1;
Alambda1_h = p.statics.Alambda1_h ;

evalLambda1 = lambda1_h(x,y,curElem,lvl,p);
evalAlambda1 = Alambda1_h(x,y,curElem,lvl,p);

val = evalAlambda1 - evalLambda1;

%% supply fourth integrand: div(A\lambda1-lambda1)(u)
function val = integrand4(x,y,curElem,lvl,p)

divLambda1_h = p.statics.divLambda1_h;
divAlambda1_h = p.statics.divAlambda1_h;
u_h = p.statics.u_h;


evalDivLambda1 = divLambda1_h(x,y,curElem,lvl,p);
evalDivAlambda1 = divAlambda1_h(x,y,curElem,lvl,p);
evalU = u_h(x,y,curElem,lvl,p);

val = (evalDivAlambda1 - evalDivLambda1).*evalU;
val = reshape(val,[1,1,length(x)]);

%% supply fith integrand: div(A\lambda1-lambda1)(u)
function val = integrand5(x,y,curElem,lvl,p)

lambda2_h = p.statics.lambda2;
Alambda2_h = p.statics.Alambda2_h;
divP_h = p.statics.divP_h;


evalLambda2 = lambda2_h(x,y,curElem,lvl,p);
evalAlambda2 = Alambda2_h(x,y,curElem,lvl,p);
evalDivP = divP_h(x,y,curElem,lvl,p);

val = (evalAlambda2 - evalLambda2).*evalDivP;
val = reshape(val,[1,1,length(x)]);

%% supply sixth integrand: div(A\p-p)(\lambda2)
function val = integrand6(x,y,curElem,lvl,p)

divP_h = p.statics.divP_h;
divAp_h = p.statics.divAp_h;
lambda2 = p.statics.lambda2;


evalDivP_h = divP_h(x,y,curElem,lvl,p);
evalDivAp_h = divAp_h(x,y,curElem,lvl,p);
evalLambda2 = lambda2(x,y,curElem,lvl,p);

val = (evalDivAp_h - evalDivP_h).*evalLambda2;
val = reshape(val,[1,1,length(x)]);

%% supply seventh integrand: (A\lambda2-lambda2)*f
function val = integrand7(x,y,curElem,lvl,p)

lambda2_h = p.statics.lambda2;
Alambda2_h = p.statics.Alambda2_h;
f = p.problem.f;

evalLambda2 = lambda2_h(x,y,curElem,lvl,p);
evalAlambda2 = Alambda2_h(x,y,curElem,lvl,p);
evalF = f(x,y,curElem,lvl,p);

val = (evalAlambda2 - evalLambda2).*evalF;
val = reshape(val,[1,1,length(x)]);
