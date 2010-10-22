function val = funcHandleStima(x,y,curElem,curLvl,p)

% author: Joscha Gedicke

kappa = p.problem.kappa;
grad  = p.statics.basisGradU;

curKappa = kappa(x,y,p);
curGrad = grad(x,y,curElem,curLvl,p);

nrBasisFuncU = size(curGrad,1);

val = zeros(nrBasisFuncU,nrBasisFuncU,length(x));

for i = 1 : length(x)    
    val(:,:,i) = curGrad(:,:,i)*curKappa(:,:,i)'*curGrad(:,:,i)'; 
end

