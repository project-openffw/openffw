function p = P1getJacobian(p)
%author: David Guenther

%% INPUT 
% load enumerated data
n4e = p.level(end).geom.n4e;
freeNodes = p.level(end).enum.freeNodes;
dofU4e = p.level(end).enum.dofU4e;
lvl = size(p.level,2);
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute the jacobian DE(x)
ST = integrate(n4e,lvl,degree,@integrand,p);
ST = permute(ST,[2 3 1]);
[I,J] = localDoFtoGlobalDoF(dofU4e,dofU4e);
jacobi = sparse(I,J,ST(:));

%% OUTPUT
p.level(end).jacobi = jacobi(freeNodes,freeNodes);

%% supply integrand: D2W(\nabla u_h)*\nabla w_h\nabla q_h
function val = integrand(x,y,curElem,lvl,p)
%% D2W(|X|)[Y,Z] = W''(|X|)*(Y*X)(Z*X) + (W'(|X|)/|X|)*(Y*Z)

DW = p.problem.nonLinearExactDer;
D2W = p.problem.nonLinearExactSecDer;
grad_h = p.statics.grad_h;
stressBasis = p.statics.stressBasis;

evalGrad = grad_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

absGrad = ( evalGrad(:,1).^2 + evalGrad(:,2).^2 ).^(1/2);
evalDW = DW(absGrad,curElem,lvl,p);
evalD2W = D2W(absGrad,curElem,lvl,p);

evalGrad = reshape(evalGrad',[1 2 length(x)]);
YZ = matMul(evalBasis,permute(evalBasis,[2 1 3]));
XY = matMul(evalGrad,permute(evalBasis,[2 1 3]));
XZ = matMul(evalGrad,permute(evalBasis,[2 1 3]));
XYXZ = matMul(permute(XY,[2 1 3]),XZ);


if norm(absGrad) > 0
    term1 = matMul(reshape(evalD2W,[1 1 length(x)]),XYXZ);
else
    term1 = 0;
end

term2 = matMul(reshape(evalDW,[1 1 length(x)]),YZ); %calculates elementwise matrix product

val = term1 + term2;
