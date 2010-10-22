function p = TWRTDWRgetJacobian(p)
% author: David Guenther, Lena Noack

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

%% compute the Jacobian DE(x) for dual problem

BT = zeros(3,3,nrElems);
CT = zeros(1,3,nrElems);

for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curEdgeLengths = length4ed(curEdges);
    
    signum = ones(1,3);
	signum(e4ed(curEdges,2) == curElem) = -1;    
    %divergence of basis functions on T,e.g. div(q) = sign*|E|/|T|
    div_qh = signum.*curEdgeLengths';
    
    %calc int_T  D^2W(Psigma_h)*sigma_h*q_h
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

%% supply integrand D^2W*(Psigma_h)*w_h*q_h
function val = integrand(x,y,curElem,lvl,p)

D2W = p.problem.conjNonLinearFuncSecDer;
Psigma_h = p.statics.sigma_h;%primal solution X
stressBasis = p.statics.stressBasis; % basis   Z

evalPSigma = Psigma_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

evalD2W = D2W(evalPSigma(:,1),evalPSigma(:,2),curElem,lvl,p);

WQ = matMul(evalBasis,permute(evalBasis,[2 1 3]));
val = matMul(reshape(evalD2W,[1 1 length(x)]),WQ);