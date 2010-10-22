function p = TWDWRgetJacobian(p)
% author: Lena Noack
%% INPUT 
n4e = p.level(end).geom.n4e;

length4ed = p.level(end).enum.length4ed;
ed4e = p.level(end).enum.ed4e;
e4ed = p.level(end).enum.e4ed;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;

lvl = size(p.level,2);

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute the Jacobian DE(x)
BT = zeros(3,3,nrElems);
CT = zeros(1,3,nrElems);
DT = zeros(3,3,nrElems);

for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curEdgeLengths = length4ed(curEdges);
    
    signum = ones(1,3);
	I = find(e4ed(curEdges,2) == curElem);
	signum(I) = -1;    
    %divergence of basis functions on T,e.g. div(q) = sign*|E|/|T|
    div_qh = signum.*curEdgeLengths';
    
    D2W = integrate(n4e(curElem,:),lvl,degree,@integrand1,p);
%     D2W = intNonLinear(@getDer4Functional,sigma_h,stressBasis,stressBasis,[],...
%                             conjNonLinearFuncDer,conjNonLinearFuncSecDer,curElem,lvl,p);
    D3W = integrate(n4e(curElem,:),lvl,degree,@integrand2,p);
%     D3W = intNonLinear(@getSecDer4Functional,sigma_h,lambda1,stressBasis,stressBasis,...
%                             conjNonLinearFuncDer,conjNonLinearFuncSecDer,curElem,lvl,p);
    %calc int_T  D^2W(p_h)*q_h*r_h
    BT(:,:,curElem) = D2W;
   	%calc int_T v_h*div(q_h)
    CT(:,:,curElem) = div_qh;	
    
    %calc int_T  D^3W(p_h)*q_h*\lambda*r_h
    DT(:,:,curElem) = D3W;
end

[I,J] = localDoFtoGlobalDoF(ed4e,ed4e);
B = sparse(I,J,BT(:));
D = sparse(I,J,DT(:));

[I,J] = localDoFtoGlobalDoF(ed4e,(1:nrElems)');
C = sparse(I,J,CT(:),nrEdges,nrElems);

dummyZeros1 = sparse(nrEdges,nrElems);
dummyZeros2 = sparse(nrElems,nrElems);
dummyZeros3 = sparse(nrEdges,nrEdges);

jacobi = [  B-D,              dummyZeros1,     -B,         -C;
            dummyZeros1',  dummyZeros2,    C',      dummyZeros2;
            B          , C            ,  dummyZeros3,  dummyZeros1;
            C',  dummyZeros2, dummyZeros1', dummyZeros2];

%% OUTPUT 
p.level(end).jacobi = jacobi;

%% supply integrand
function val = integrand1(x,y,curElem,lvl,p)
F1 = p.problem.F1;
F2 = p.problem.F2;
CONV = p.params.CONV;

if strcmp(CONV,'c')
%% D2W**(F)[G,H] = W''**1(F)(G,H) + W''**2(F)(F,G)(F,H) + W''**3(F)(F2,H)(F2,G)
D2W1 = p.problem.nonLinearRegSecDerA;
D2W2 = p.problem.nonLinearRegSecDerB;
D2W3 = p.problem.nonLinearRegSecDerC;
else
% D2W(F)[G,H] = W''1(F)(G,H) + W''2(F)((F-F1,G)*(F-F2,H) + (F-F1,H)*(F-F2,G))
D2W1 = p.problem.nonLinearExactSecDerA;
D2W2 = p.problem.nonLinearExactSecDerB;
end

grad_h = p.statics.grad_h;
stressBasis = p.statics.stressBasis;
basisU = p.statics.basisU;

evalGrad = grad_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);
evalBasisU = basisU(x,y,curElem,lvl,p);

evalBasisU=reshape(evalBasisU,[3 1 length(x)]);
evalBasisU2=permute(evalBasisU,[2 1 3]);
evalBasisU2=reshape(evalBasisU2,[1 3 length(x)]);

wq = matMul(evalBasisU,evalBasisU2);

grad1 = evalGrad(:,1);
grad2 = evalGrad(:,2);
evalGrad = reshape(evalGrad',[1 2 length(x)]);

evalD2W1 = D2W1(grad1,grad2,curElem,lvl,p);
evalD2W2 = D2W2(grad1,grad2,curElem,lvl,p);
if strcmp(CONV,'c')
evalD2W3 = D2W3(grad1,grad2,curElem,lvl,p);
end

GH = matMul(evalBasis,permute(evalBasis,[2 1 3]));
evalGrad1 = reshape([grad1-F1(1) grad2-F1(2)]',[1 2 length(x)]);
evalGrad2 = reshape([grad1-F2(1) grad2-F2(2)]',[1 2 length(x)]);
F2vec = reshape([F2(1)*ones(size(grad1,1),1) F2(2)*ones(size(grad1,1),1)]',[1 2 length(x)]);
F1G = matMul(evalGrad1,permute(evalBasis,[2 1 3]));
F2G = matMul(evalGrad2,permute(evalBasis,[2 1 3]));
F1H = matMul(evalGrad1,permute(evalBasis,[2 1 3]));
F2H = matMul(evalGrad2,permute(evalBasis,[2 1 3]));
FG = matMul(evalGrad,permute(evalBasis,[2 1 3]));
FH = matMul(evalGrad,permute(evalBasis,[2 1 3]));
GF2 = matMul(F2vec,permute(evalBasis,[2 1 3]));
HF2 = matMul(F2vec,permute(evalBasis,[2 1 3]));
F1GF2H = matMul(permute(F1G,[2 1 3]),F2H);
F1HF2G = matMul(permute(F1H,[2 1 3]),F2G);
FGFH = matMul(permute(FG,[2 1 3]),FH);
GF2HF2 = matMul(permute(GF2,[2 1 3]),HF2);

%% D2W**(F)[G,H] = W''**1(F)(G,H) + W''**2(F)(F,G)(F,H) + W''**3(F)(F2,G)(F2,H)
% D2W(F)[G,H] = W''1(F)(G,H) + W''2(F)((F-F1,G)*(F-F2,H) + (F-F1,H)*(F-F2,G))

if strcmp(CONV,'c')
term1 = matMul(reshape(evalD2W1,[1 1 length(x)]),GH)+matMul(reshape(evalD2W2,[1 1 length(x)]),FGFH);
term2 = matMul(reshape(evalD2W3,[1 1 length(x)]),GF2HF2);
else
term1 = matMul(reshape(evalD2W1,[1 1 length(x)]),GH);
term2 = matMul(reshape(evalD2W2,[1 1 length(x)]),F1GF2H+F1HF2G);
end

val = term1 + term2 + 2*wq;

%% supply integrand D3W
function val = integrand2(x,y,curElem,lvl,p)
% computation of D3W(\nabla u_h)*\nabla w_h\nabla q_h\nabla p_h
% D3W**(F) = W''**2(F)[(F,J)(G,H) + (F,H)(G,J) + (F,G)(H,J)]
% D3W  (F) = 4(G,H)(J,2F-F1-F2) + 4(G,J)(H,2F-F1-F2) + 4(H,J)(G,2F-F1-F2)

F1 = p.problem.F1;
F2 = p.problem.F2;
CONV = p.params.CONV;

if strcmp(CONV,'c')
%% D3W**(F) = W''**2(F)[(F,J)(G,H) + (F,H)(G,J) + (F,G)(H,J)]
D2W = p.problem.nonLinearRegSecDerB;
else
% D3W  (F) = 4(G,H)(J,2F-F1-F2) + 4(G,J)(H,2F-F1-F2) + 4(H,J)(G,2F-F1-F2)
end

grad_h = p.statics.grad_h;
evalGrad = grad_h(x,y,curElem,lvl,p);
grad1 = evalGrad(:,1);
grad2 = evalGrad(:,2);
evalGrad = reshape(evalGrad',[1 2 length(x)]);


sigma_h = p.statics.sigma_h;
lambda1 = p.statics.lambda1;
stressBasis = p.statics.stressBasis;

evalSigma = sigma_h(x,y,curElem,lvl,p);
evalLambda = lambda1(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

absSigma = ( evalSigma(:,1).^2 + evalSigma(:,2).^2 ).^(1/2);

if strcmp(CONV,'c')
evalDW = D2W(absSigma,curElem,lvl,p);
end

evalSigma = reshape(evalSigma',[1 2 length(x)]);
evalLambda = reshape(evalLambda',[1 2 length(x)]);
FG = matMul(evalSigma,permute(evalBasis,[2 1 3]));
FH = matMul(evalSigma,permute(evalBasis,[2 1 3]));
FJ = matMul(evalSigma,permute(evalLambda,[2 1 3]));
JG = matMul(evalLambda,permute(evalBasis,[2 1 3]));
JH = matMul(evalLambda,permute(evalBasis,[2 1 3]));
GH = matMul(evalBasis,permute(evalBasis,[2 1 3]));
FJGH = matMul(permute(FJ,[2 1 3]),GH);
FHGJ = matMul(permute(FH,[2 1 3]),JG);
FGHJ = matMul(permute(FG,[2 1 3]),JH);


evalGrad = reshape([2*grad1-F1(1)-F2(1) 2*grad2-F1(2)-F2(2)]',[1 2 length(x)]); %(2F-F1-F2)
diffG = matMul(evalGrad,permute(evalBasis,[2 1 3]));
diffH = matMul(evalGrad,permute(evalBasis,[2 1 3]));
diffJ = matMul(evalGrad,permute(evalLambda,[2 1 3]));

dGHJ = matMul(permute(diffG,[2 1 3]),JH);
dHGJ = matMul(permute(diffH,[2 1 3]),JG);
dJGH = matMul(permute(diffJ,[2 1 3]),GH);

% D3W**(F) = W''**2(F)[(F,J)(G,H) + (F,H)(G,J) + (F,G)(H,J)]
% D3W  (F) = 4(G,H)(J,2F-F1-F2) + 4(G,J)(H,2F-F1-F2) + 4(H,J)(G,2F-F1-F2)

if strcmp(CONV,'c')
val = matMul(reshape(evalDW,[1 1 length(x)]),FJGH+FHGJ+FGHJ);
else
val = 4*(dGHJ+dHGJ+dJGH);
end
