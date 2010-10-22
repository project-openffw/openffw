function val = funcHandleRHSVolume(x,y,curElem,curLvl,p)
% function handle for f*phi_j
% use: val = funcHandleRHSVolume(x,y,curElem,curlvl,p)

%% INPUT
f = p.problem.f;
basisU = p.statics.basisU;
dofU4e = p.level(curLvl).enum.dofU4e;
nrBasisFuncU = size(dofU4e,2);

%% f*phi_j
evalF = f(x,y,curElem,curLvl,p);
dimf = size(evalF,2);
evalBasisU = basisU(x,y,curElem,curLvl,p);
evalBasisU = reshape(evalBasisU, [length(x) dimf nrBasisFuncU]);
evalF = evalF(:)*ones(1,nrBasisFuncU);
evalF = reshape(evalF,[length(x) dimf nrBasisFuncU ]);
integrand = evalF .* evalBasisU;
integrand = sum(integrand,2);

val = permute(integrand, [ 3 2 1 ]);
