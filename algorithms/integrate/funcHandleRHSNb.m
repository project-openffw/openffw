function val = funcHandleRHSNb(x,y,curEd,curLvl,p)
% function handle for f*phi_j
% use: val = funcHandleRHSNb(x,y,curElem,curlvl,p)

%% INPUT
g = p.problem.g;
basisU = p.statics.basisU;
dofU4e = p.level(curLvl).enum.dofU4e;
nrBasisFuncU = size(dofU4e,2);
e4ed = p.level(curLvl).enum.e4ed;
NbEd = p.level(curLvl).enum.NbEd;

c4n = p.level(curLvl).geom.c4n;
Nb  = p.level(curLvl).geom.Nb;
index = find(NbEd == curEd);
normal = getNormals4NbEd(c4n,Nb(index,:));

%% f*phi_j
evalG = g(x,y,(normal(:)*ones(1,length(x)))',p);
dimg = size(evalG,2);
evalBasisU = basisU(x,y,e4ed(curEd),curLvl,p);
evalBasisU = reshape(evalBasisU, [length(x) dimg nrBasisFuncU]);
evalG = evalG(:)*ones(1,nrBasisFuncU);
evalG = reshape(evalG,[length(x) dimg nrBasisFuncU ]);
integrand = evalG .* evalBasisU;
integrand = sum(integrand,2);

val = permute(integrand, [ 3 2 1 ]);
