function p = TWP2getFuncVal(p)
% author: David Guenther, Lena Noack

%% INPUT
% load enumerated data
n4e = p.level(end).geom.n4e;
dofU4e = p.level(end).enum.dofU4e;
f4e = p.level(end).f4e;
lvl = size(p.level,2);
freeNodes = p.level(end).enum.freeNodes;
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute the function value E(x)
ST = integrate(n4e,lvl,degree,@integrand,p);

I = dofU4e;
S = accumarray(I(:),ST(:));
rhs = accumarray(I(:),f4e(:));

funcVal = S - rhs;

%% OUTPUT 
p.level(end).funcVal = funcVal(freeNodes);

%% supply integrand: DW(\nabla u_h)*\nabla w_h + 2(u_h - f)w_h
function val = integrand(x,y,curElem,lvl,p)

u_h = p.statics.u_h;
f = p.problem.f0;
evalUh = u_h(x,y,curElem,lvl,p);
evalF = f(x,y,curElem,lvl,p);
evalDif = reshape(2*(evalUh-evalF)',[1 1 length(x)]);

sigma_h = p.statics.sigma_h;
stressBasis = p.statics.stressBasis;

evalSigma = sigma_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

evalSigma = reshape(evalSigma',[2 1 length(x)]);

rhs = matMul(evalDif,evalBasis); %rhs = 2(u_h-f)w_h

val = matMul(evalBasis,evalSigma) + rhs;
