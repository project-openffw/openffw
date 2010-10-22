function p = TWRTgetJacobian(p)
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

x = p.level(end).x0;
lvl = size(p.level,2);
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute the Jacobian DE(x)

BT = zeros(3,3,nrElems);
CT = zeros(1,3,nrElems);
DT = zeros(1,1,nrElems);

for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curEdgeLengths = length4ed(curEdges);
    curX4ed = x(curEdges);
    curX4e = x(nrEdges + curElem);
    
    signum = ones(1,3);
	signum(e4ed(curEdges,2) == curElem) = -1;    
    %divergence of basis functions on T,e.g. div(q) = sign*|E|/|T|
    div_qh = signum.*curEdgeLengths';
    
    %calc int_T  D^2W(p_h)*w_h*q_h
    D2W = integrate(n4e(curElem,:),lvl,degree,@integrand,p);                   
    BT(:,:,curElem) = D2W;
   	%calc int_T v_h*div(q_h)
    CT(:,:,curElem) = div_qh;	
   	%calc -2 int_T v_h*v_h
    DT(:,:,curElem) = -2* curX4ed'*curX4ed;	
end

[I,J] = localDoFtoGlobalDoF(dofSigma4e,dofSigma4e);
B = sparse(I,J,BT(:));

[I,J] = localDoFtoGlobalDoF(dofSigma4e,dofU4e);
C = sparse(I,J,CT(:),nrEdges,nrElems);

[I,J] = localDoFtoGlobalDoF(dofU4e,dofU4e);
D = sparse(I,J,DT(:));

jacobi = [  B,    C
%            C',   sparse(length(dofU4e),length(dofU4e))];
            C',   D ];

%% OUTPUT 
p.level(end).jacobi = jacobi;

%% supply integrand D^2W*(sigma_h)*w_h*q_h
function val = integrand(x,y,curElem,lvl,p)
D2W = p.problem.conjNonLinearFuncSecDer;

sigma_h = p.statics.sigma_h;
stressBasis = p.statics.stressBasis;

evalSigma = sigma_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);
%WQ = matMul(evalBasis,permute(evalBasis,[2 1 3]));

evalD2W = D2W(evalSigma(:,1),evalSigma(:,2),curElem,lvl,p);

%val = matMul(reshape(evalD2W,[1 1 length(x)]),WQ);
Q_D2W = matMul(evalBasis,evalD2W);
val = matMul(Q_D2W,permute(evalBasis,[2 1 3]));
val = full(val);