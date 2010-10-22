function p = ODP1estimate_Dual(p)
%author: Lena Noack
% use Repin's duality theory

%% INPUT
n4e = p.level(end).geom.n4e;
lvl = size(p.level,2);
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);
nrElems = p.level(end).nrElems;
%dofU4e = p.level(end).enum.dofU4e;
%f4e = p.level(end).f4e;
midPoint4e = p.level(end).enum.midPoint4e;


gradh = p.statics.GradU_h;
sigmah = p.statics.sigma_h;
stressBasis = p.statics.stressBasis;
DW_eps = p.problem.conjNonLinearFuncDer;

grad_h = zeros(nrElems,2);
sigma_h = zeros(nrElems,2);
tau = zeros(nrElems,2);
for curElem = 1:nrElems
    grad_h(curElem,:) = gradh(midPoint4e(curElem,1),midPoint4e(curElem,2),curElem,lvl,p);
    sigma_h(curElem,:) = sigmah(midPoint4e(curElem,1),midPoint4e(curElem,2),curElem,lvl,p);

%    tau(curElem,:) = [midPoint4e(curElem,1) 0];
    tau(curElem,:) = 0.5*[midPoint4e(curElem,1) midPoint4e(curElem,2)];

    evalBasis = stressBasis(midPoint4e(curElem,1),midPoint4e(curElem,2),curElem,lvl,p);

    absTau = ( tau(curElem,1).^2 + tau(curElem,2).^2 ).^(1/2);

    evalDW_eps = DW_eps(absTau,curElem,lvl,p);
    evalDW(curElem,:) = sum(evalDW_eps*evalBasis,1);
end

term1 = integrate(n4e,lvl,degree,@integrandW,p);

term2 = integrate(n4e,lvl,degree,@integrandWAst,p);


term3 = matMul(reshape(sigma_h,[1 2 nrElems]),reshape(grad_h,[2 1 nrElems]));
term3 = reshape(term3,[nrElems 1]);

term4 = (1/2)*matMul(reshape(sigma_h-tau,[1 2 nrElems]),reshape(evalDW - grad_h,[2 1 nrElems]));
term4 = reshape(term4,[nrElems 1]);

eta4T = term1 + term2 + term3 + term4;
%eta4T = term1 + term2 + term3 + term4;

%% OUTPUT
p.level(end).etaT = eta4T;
p.level(end).estimatedError = abs(sum(eta4T));  %est. for ||sigma-sigma_h||

% W(Du_h)
function val = integrandW(x,y,curElem,lvl,p)
W = p.problem.nonLinearExact;
grad_h = p.statics.GradU_h;
evalGrad = grad_h(x,y,curElem,lvl,p);
absGrad = ( evalGrad(:,1).^2 + evalGrad(:,2).^2 ).^(1/2);
evalW = W(absGrad,curElem,lvl,p);
val = reshape(evalW,[1 1 length(x)]);

% W*(-sigma_h)
function val = integrandWAst(x,y,curElem,lvl,p)
sigma_h = p.statics.sigma_h;
evalSigma = sigma_h(x,y,curElem,lvl,p);
absSigma = -( evalSigma(:,1).^2 + evalSigma(:,2).^2 ).^(1/2);

W_eps = p.problem.conjNonLinearFunc;
evalW = W_eps(absSigma,curElem,lvl,p);
val = reshape(evalW,[1 1 length(x)]);

