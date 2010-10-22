function p = ODRTgetJacobian(p)
% author: David Guenther 

%% INPUT 
n4e = p.level(end).geom.n4e;
length4ed = p.level(end).enum.length4ed;
ed4e = p.level(end).enum.ed4e;
e4ed = p.level(end).enum.e4ed;
dofU4e = p.level(end).enum.dofU4e;
dofSigma4e = p.level(end).enum.dofSigma4e;

nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;

lvl = size(p.level,2);
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute the Jacobian DE(x)

BT = zeros(3,3,nrElems);
CT = zeros(1,3,nrElems);

for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curEdgeLengths = length4ed(curEdges);
    
    signum = ones(1,3);
	signum(e4ed(curEdges,2) == curElem) = -1;    
    %divergence of basis functions on T,e.g. div(q) = sign*|E|/|T|
    div_qh = signum.*curEdgeLengths';
    
    %calc int_T  D^2W(p_h)*w_h*q_h
    D2W = integrate(n4e(curElem,:),lvl,degree,@integrand,p);                   
    BT(:,:,curElem) = D2W;
   	%calc int_T v_h*div(q_h)
    CT(:,:,curElem) = div_qh;	
end

[I,J] = localDoFtoGlobalDoF(dofSigma4e,dofSigma4e);
B = sparse(I,J,BT(:));

[I,J] = localDoFtoGlobalDoF(dofSigma4e,dofU4e);
C = sparse(I,J,CT(:),nrEdges,nrElems);

jacobi = [  B,    C
            C',  sparse(length(dofU4e),length(dofU4e))];

%% OUTPUT 
p.level(end).jacobi = jacobi;

%% supply integrand
function val = integrand(x,y,curElem,lvl,p)
% W''(|X|)/|X|^2*(X*Y)(X*Z) + W'(|X|)/|X|^3*( Y*Z*|X|^2 - (X*Y)(X*Z) )

DW = p.problem.conjNonLinearFuncDer;
D2W = p.problem.conjNonLinearFuncSecDer;
sigma_h = p.statics.sigma_h;
stressBasis = p.statics.stressBasis;

evalSigma = sigma_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

absSigma = ( evalSigma(:,1).^2 + evalSigma(:,2).^2 ).^(1/2);
evalDW = DW(absSigma,curElem,lvl,p);
evalD2W = D2W(absSigma,curElem,lvl,p);

evalSigma = reshape(evalSigma',[1 2 length(x)]);
YZ = matMul(evalBasis,permute(evalBasis,[2 1 3]));
XY = matMul(evalSigma,permute(evalBasis,[2 1 3]));
XZ = matMul(evalSigma,permute(evalBasis,[2 1 3]));
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