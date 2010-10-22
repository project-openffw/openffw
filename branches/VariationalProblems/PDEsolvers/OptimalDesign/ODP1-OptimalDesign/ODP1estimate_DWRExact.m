function p = ODP1estimate_DWR(p)
%author: Lena Noack

%% INPUT
n4e = p.level(end).geom.n4e;
lvl = size(p.level,2);
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);
nrElems = p.level(end).nrElems;
%dofU4e = p.level(end).enum.dofU4e;
%f4e = p.level(end).f4e;
DWRbool = loadField('p.params','DWRbool',p,0);
if DWRbool == 1
      GOAL = p.params.GOAL;

      %GOAL_IM: mean integral:   J(tau,v) = 1/|Omega1|  int_Omega1 tau dx
      if strcmp(GOAL,'GOAL_IM')
        GoalFunc = integrate(n4e,lvl,degree,@integrandJMid,p);

      %GOAL_DP: dual problem:      J(tau,v) = int_Omega W*_epsilon(tau) dx
      elseif strcmp(GOAL,'GOAL_DP')
        GoalFunc = integrate(n4e,lvl,degree,@integrandJDW,p);

      else
        GoalFunc = integrate(n4e,lvl,degree,@integrandJPR,p);

      end
      
      
      D2W = integrate(n4e,lvl,degree,@integrandD2W,p);
%      Alpha = integrate(n4e,lvl,degree,@integrandAlpha,p); %only for 2-Well Benchmark
      DW = integrate(n4e,lvl,degree,@integrandDW,p);
%      Alpha2 = integrate(n4e,lvl,degree,@integrandAlpha2,p); %only for 2-Well Benchmark
      RHS = integrate(n4e,lvl,degree,@integrandF,p);

%      eta4T = GoalFunc - D2W + Alpha - DW - Alpha2 + RHS; %only for 2-Well Benchmark
      %eta4TNodes = GoalFunc - D2W - DW + RHS;
      %eta4T = sum(eta4TNodes,1);
%      eta4T = reshape(sqrt(matMul(reshape(eta4TNodes,[1 3 nrElems]),reshape(eta4TNodes,[3 1 nrElems]))),[1 nrElems]);
      rhoD = GoalFunc - D2W;
      rhoP = -DW + RHS;
else
    fprintf('No dual estimator calculated, select goal function and set DWRbool=1.');
    eta4T = zeros(nrElems,1);
end

%eta4T=0.5*eta4T;
residual = (1/2*rhoD + 1/2*rhoP).^2;
%residual = (1/2*rhoD - 1/2*rhoP).^2;
%linError = (1/2*(rhoD - rhoP)).^2; %-> eta = 0.5(rho* + rho) + Delta rho
linError = 0;                       %-> eta = 0.5(rho* + rho)

etaT = sqrt(abs(residual));
etaOsc = sqrt(abs(linError));
estError = sqrt(abs(sum(residual) + sum(linError)));

%% OUTPUT
p.level(end).etaT = etaT;
p.level(end).etaOsc = etaOsc;
p.level(end).estimatedError = estError;
%p.level(end).estimatedError = sqrt(abs(sum(eta4T)));  %est. for ||sigma-sigma_h||

function val = integrandJMid(x,y,curElem,lvl,p)

%basis = p.statics.u_h;
basis = p.statics.basisU;
u_exact = p.problem.u_exact;

evalUexact = u_exact(x,y,p);
evalBasis = basis(x,y,curElem,lvl,p);

evalBasis = evalUexact*ones(1,3) - evalBasis;
nrBasis = size(evalBasis,2);

evalF = zeros(size(x,1),nrBasis);
%f=1 on Omega2 (Lshape domain) and 0 on O1,O3:
%  _________
% |    |    |
% | O1 | O2 |
% |____|____|
% |    |
% | O3 |
% |____|
%
%index1 = find( x >= 0 );
%evalF(index1,:) = 1*C;

%f=1 on Omega2 (SquareSlit domain) and 0 on O1:
%  _____
% |  |  |
% |O1|O2|
% |__|__|
%
index1 = find( x >= 0 );
evalF(index1,:) = 0.5;

val = reshape(evalF'.*evalBasis',[nrBasis 1 length(x)]);
val = sum(val,1);

% supply integrand int_O1 DW*_varepsilon(sigma_h,tau_h)
function val = integrandJDW(x,y,curElem,lvl,p)

grad_h = p.statics.grad_h;
basis = p.statics.stressBasis;
u_exact = p.problem.u_exact;

evalUexact = u_exact(x,y,p);
evalGrad = grad_h(x,y,curElem,lvl,p);
evalBasis = basis(x,y,curElem,lvl,p);
evalUexact = matMul(reshape(ones(6,length(x)),[3 2 length(x)]),reshape(evalUexact,[1 1 length(x)]));

evalBasis = evalUexact - evalBasis;
nrBasis = size(evalBasis,2);

evalF = zeros(size(x,1),3);
index1 = find( x <= -0.5 );
evalF(index1,:) = 1;

evalGrad = reshape(evalGrad,[1 2 length(x)]);   
evalF = reshape(evalF,[1 3 length(x)]);   
val = matMul(evalGrad,permute(evalBasis,[2 1 3])); 
val = matMul(evalF,permute(val,[2 1 3])); 


% supply int_Omega DW(Du_h) D2W(Du_h,Du-Dv_h) dx
function val = integrandJPR(x,y,curElem,lvl,p)

DW = p.problem.nonLinearExactDer;
D2W = p.problem.nonLinearExactSecDer;
grad_h = p.statics.grad_h;
stressBasis = p.statics.stressBasis;
Du_exact = p.problem.gradU_exact;

evalGrad = grad_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);
evalDUexact = Du_exact(x,y,p);

evalDUexact = reshape(evalDUexact,[1,2,length(x)]);
evalDU_exact(1,:,:) = evalDUexact;
evalDU_exact(2,:,:) = evalDUexact;
evalDU_exact(3,:,:) = evalDUexact;
evalBasis = evalDU_exact-evalBasis;

absGrad = ( evalGrad(:,1).^2 + evalGrad(:,2).^2 ).^(1/2);
evalDW = DW(absGrad,curElem,lvl,p);                     % D'W(Du)/Du
evalD2W = D2W(absGrad,curElem,lvl,p);                   % D''W(Du)


evalBasis = reshape(evalBasis,[3 2 length(x)]);
evalDW = reshape(evalDW,[1 1 length(x)]);            % 1 1 x           
evalD2W = reshape(evalD2W,[1 1 length(x)]);          % 1 1 x

if norm(absGrad) > 0
    evalSigGrad(:,1) = evalGrad(:,1)./absGrad;
    evalSigGrad(:,2) = evalGrad(:,2)./absGrad;
else
    evalSigGrad = zeros(2,length(x));
end
evalSigGrad = reshape(evalSigGrad,[1 2 length(x)]);
XY = matMul(evalSigGrad,permute(evalBasis,[2 1 3])); % 1 3 x
XYX = matMul(permute(XY,[2 1 3]),evalSigGrad);       % 3 2 x

D2W_Vh = matMul(evalD2W,XYX);                        % 3 2 x
DW_Vh  = matMul(evalDW,evalBasis - XYX);             % 3 2 x

evalGrad = reshape(evalGrad,[1 2 length(x)]);
term1 = matMul(evalDW,evalGrad);                     % 1 2 x
term2 = D2W_Vh + DW_Vh;                              % 3 2 x

val = 2*matMul(term1,permute(term2,[2 1 3]));        % 1 3 x
val = permute(val,[2 1 3]);

val = sum(val,1);

%% supply integrand: D2W(Du_h)*(Du-Dv_h)Dz_h
function val = integrandD2W(x,y,curElem,lvl,p)

% W''(|X|)/|X|^2*(X*Y)(X*Z) + W'(|X|)/|X|^3*( Y*Z*|X|^2 - (X*Y)(X*Z) )

DW = p.problem.nonLinearExactDer;
D2W = p.problem.nonLinearExactSecDer;
grad_h = p.statics.grad_h;
DWRgrad_h = p.statics.DWRgrad_h;
stressBasis = p.statics.stressBasis;
Du_exact = p.problem.gradU_exact;

evalDUexact = Du_exact(x,y,p);
evalGrad = grad_h(x,y,curElem,lvl,p);
evalDWRGrad = DWRgrad_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

evalDUexact = reshape(evalDUexact,[1,2,length(x)]);
evalDU_exact(1,:,:) = evalDUexact;
evalDU_exact(2,:,:) = evalDUexact;
evalDU_exact(3,:,:) = evalDUexact;
evalBasis = evalDU_exact-evalBasis;

absGrad = ( evalGrad(:,1).^2 + evalGrad(:,2).^2 ).^(1/2);
evalDW = DW(absGrad,curElem,lvl,p);
evalD2W = D2W(absGrad,curElem,lvl,p);

% W''(|X|)/|X|^2*(X*Y)(X*Z) + W'(|X|)/|X|^3*( Y*Z*|X|^2 - (X*Y)(X*Z) )
evalGrad = reshape(evalGrad',[1 2 length(x)]);
evalDWRGrad = reshape(evalDWRGrad',[1 2 length(x)]);
YZ = matMul(evalBasis,permute(evalDWRGrad,[2 1 3]));
XY = matMul(evalGrad,permute(evalBasis,[2 1 3]));
XZ = matMul(evalGrad,permute(evalDWRGrad,[2 1 3]));
XYXZ = matMul(permute(XY,[2 1 3]),XZ);

if norm(absGrad) > 0
    term1 = matMul(reshape(evalD2W./absGrad.^2,[1 1 length(x)]),XYXZ);
    term2 = -matMul(reshape(evalDW./absGrad.^2,[1 1 length(x)]),XYXZ);
else
    term1 = 0;
    term2 = 0;
end

term3 = matMul(reshape(evalDW,[1 1 length(x)]),YZ); %calculates elementwise matrix product

val = term1 + term2 + term3;
val = reshape(val,[3 1 length(x)]);
val = sum(val,1);

%% supply integrand: DW(Du_h)*(Dz-Dv_h)
function val = integrandDW(x,y,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
stressBasis = p.statics.stressBasis;
DWRgrad_h = p.statics.DWRgrad_h;

evalSigma = sigma_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);
evalDWRGrad = DWRgrad_h(x,y,curElem,lvl,p);

evalDWRGrad = reshape(evalDWRGrad,[1,2,length(x)]);
evalDWR_Grad(1,:,:) = evalDWRGrad;
evalDWR_Grad(2,:,:) = evalDWRGrad;
evalDWR_Grad(3,:,:) = evalDWRGrad;
evalBasis = evalDWR_Grad-evalBasis;

evalSigma = reshape(evalSigma,[2 1 length(x)]);

val = matMul(evalBasis,evalSigma); 
val = sum(val,1);

%% supply integrand: f(z-v_h)
function val = integrandF(x,y,curElem,lvl,p)
f = p.problem.f;
basis = p.statics.basisU;
DWRu_h = p.statics.DWRu_h;

evalBasis = basis(x,y,curElem,lvl,p);
evalDWRu_h = DWRu_h(x,y,curElem,lvl,p);
nrBasis = size(evalBasis,2);

evalF = f(x,y,curElem,lvl,p)*ones(1,nrBasis);

evalDWRuh(:,1) = evalDWRu_h;
evalDWRuh(:,2) = evalDWRu_h;
evalDWRuh(:,3) = evalDWRu_h;

val = reshape(evalF'.*(evalDWRuh-evalBasis)',[nrBasis 1 length(x)]);
val = sum(val,1);
