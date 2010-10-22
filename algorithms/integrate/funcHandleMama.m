function val = funcHandleMama(x,y,curElem,curLvl,p)

% author: Joscha Gedicke

mu = p.problem.mu;
basis = p.statics.basisU;

curMu = mu(x,y,p);
curBasis = basis(x,y,curElem,curLvl,p);

nrBasisFuncU = size(curBasis,2);

val = zeros(nrBasisFuncU,nrBasisFuncU,length(x));

for i = 1 : length(x)    
    val(:,:,i) = curMu(i)*curBasis(i,:)'*curBasis(i,:);    
end

