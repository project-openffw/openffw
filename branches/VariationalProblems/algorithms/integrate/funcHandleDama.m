function val = funcHandleDama(x,y,curElem,curLvl,p)

% author: Joscha Gedicke

lambda = p.problem.lambda;
grad  = p.statics.basisGradU;
basis = p.statics.basisU;

curLambda = lambda(x,y,p);
curGrad = grad(x,y,curElem,curLvl,p);
curBasis = basis(x,y,curElem,curLvl,p);

nrBasisFuncU = size(curGrad,1);

val = zeros(nrBasisFuncU,nrBasisFuncU,length(x));

for i = 1 : length(x)    
    val(:,:,i) = curGrad(:,:,i)*curLambda(i,:)'*curBasis(i,:);    
end

