function p = ODRTDWRgetFuncVal(p)
% author: David Guenther, Lena Noack

%% INPUT
n4e = p.level(end).geom.n4e;
length4ed = p.level(end).enum.length4ed;
ed4e = p.level(end).enum.ed4e;
e4ed = p.level(end).enum.e4ed;
n4ed = p.level(end).enum.n4ed;
DbEd = p.level(end).enum.DbEd;
dofU4e = p.level(end).enum.dofU4e;
dofSigma4e = p.level(end).enum.ed4e;%p.level(end).enum.dofSigma4e;

nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;
nrNodes = p.level(end).nrNodes;


x = p.level(end).DWRx0;
lvl = size(p.level,2);

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute the function value E(x) for dual problem

BT = zeros(1,3,nrElems);
CT = zeros(1,3,nrElems);
DT = zeros(1,1,nrElems);
FT = zeros(1,3,nrElems);

for curElem = 1:nrElems    
	curEdges = ed4e(curElem,:);
	curEdgeLengths = length4ed(curEdges);
    curX4ed = x(curEdges);
    curX4e = x(nrEdges + curElem);
    
    signum = ones(1,3);
	signum(e4ed(curEdges,2) == curElem) = -1;    
    %divergence of basis functions on T,e.g. div(q) = sign*|E|/|T|
    div_qh = signum.*curEdgeLengths';
    
    %calc int_T  D2W*(Psigma_h)*sigma_h*q_h
    D2W = integrate(n4e(curElem,:),lvl,degree,@integrand,p);
    F_qh = integrate(n4e(curElem,:),lvl,degree,@integrandF,p);
    
    BT(:,:,curElem) = D2W;
   	%calc int_T u_h*div(q_h)
    CT(:,:,curElem) = curX4e*div_qh;	
    %calc int_T v_h*div(p_h)
    DT(:,:,curElem) = curX4ed'*div_qh';	
    FT(:,:,curElem) = -F_qh;
end

I = dofSigma4e';
B = accumarray(I(:),BT(:));
C = accumarray(I(:),CT(:));
F = accumarray(I(:),FT(:));
D = accumarray(dofU4e(:),DT(:));

boundary = zeros(nrEdges + nrElems,1);
boundary(DbEd) = integrate(n4ed(DbEd,:),lvl,degree,@integrandBoundary,p);

b = -boundary;
f4e = p.level(end).DWRf4e;
b(nrEdges+1:end) = f4e;

funcVal = [B + C + F;
               D]   +   b;

%% OUTPUT 
p.level(end).funcVal = funcVal;

%% supply integrand D^2W*(Psigma_h;Dsigma_h;phi_h)
function val = integrand(x,y,curElem,lvl,p)
% W''(|X|)/|X|^2*(X*Y)(X*Z) + W'(|X|)/|X|^3*( Y*Z*|X|^2 - (X*Y)(X*Z) )

DW = p.problem.conjNonLinearFuncDer;
D2W = p.problem.conjNonLinearFuncSecDer;
Psigma_h = p.statics.sigma_h;     %primal solution X
Dsigma_h = p.statics.DWRsigma_h;  %dual solution   Y
stressBasis = p.statics.stressBasis;     % basis   Z

evalPSigma = Psigma_h(x,y,curElem,lvl,p);
evalDSigma = Dsigma_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

absPSigma = ( evalPSigma(:,1).^2 + evalPSigma(:,2).^2 ).^(1/2);
evalDW = DW(absPSigma,curElem,lvl,p);
evalD2W = D2W(absPSigma,curElem,lvl,p);

evalPSigma = reshape(evalPSigma',[1 2 length(x)]);
evalDSigma = reshape(evalDSigma',[1 2 length(x)]);
YZ = matMul(evalDSigma,permute(evalBasis,[2 1 3]));
XY = matMul(evalPSigma,permute(evalDSigma,[2 1 3]));
XZ = matMul(evalPSigma,permute(evalBasis,[2 1 3]));
XYXZ = matMul(permute(XY,[2 1 3]),XZ);

if norm(absPSigma) > 0
    term1 = matMul(reshape(evalD2W./absPSigma.^2,[1 1 length(x)]),XYXZ);
    term2 = -matMul(reshape(evalDW./absPSigma.^2,[1 1 length(x)]),XYXZ);
else
    term1 = 0;
    term2 = 0;
end

term3 = matMul(reshape(evalDW,[1 1 length(x)]),YZ);

val = term1 + term2 + term3;

%% supply integrand for Dirichlet boundary
function val = integrandBoundary(x,y,curEdge,lvl,p)

u_D = p.problem.u_D;
stressBasis = p.statics.stressBasis;
e4ed = p.level(lvl).enum.e4ed;
ed4e = p.level(lvl).enum.ed4e;
normals4ed = p.level(lvl).enum.normals4ed;

curElem = e4ed(curEdge,1);
edges = ed4e(curElem,:);
index = find(edges == curEdge);
normal = normals4ed(curEdge,:);

evalU_D = u_D(x,y,p);

evalStressBasis = stressBasis(x,y,curElem,lvl,p);
curBasisFunc = squeeze(evalStressBasis(index,:,:))';
basis4normal = curBasisFunc*normal';

val = evalU_D.*basis4normal;
val = reshape(val,[1 1 length(x)]);

%% supply f_DWR * q_h (first goal function)
function val = integrandF(x,y,curElem,lvl,p)

%stressBasis = p.statics.stressBasis; 
%evalBasis = stressBasis(x,y,curElem,lvl,p);
%
%f4e = p.level(end).DWRf4e;
%evalF4e = f4e(curElem);
%F4e = evalF4e*ones(2,length(x));
%F4e = reshape(F4e,[1 2 length(x)]);
%
%FQ = matMul(F4e,permute(evalBasis,[2 1 3]));
%
%val = FQ;

DW = p.problem.nonLinearExactDer;
D2W = p.problem.nonLinearExactSecDer;
grad_h = p.statics.grad_h;
DWRgrad_h = p.statics.DWRgrad_h;
stressBasis = p.statics.stressBasis;

evalGrad = grad_h(x,y,curElem,lvl,p);
evalDWRGrad = DWRgrad_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

absGrad = ( evalGrad(:,1).^2 + evalGrad(:,2).^2 ).^(1/2);
evalDW = DW(absGrad,curElem,lvl,p);
evalD2W = D2W(absGrad,curElem,lvl,p);

evalGrad = reshape(evalGrad',[1 2 length(x)]);
evalDWRGrad = reshape(evalDWRGrad',[1 2 length(x)]);
YZ = matMul(evalDWRGrad,permute(evalBasis,[2 1 3]));
XY = matMul(evalGrad,permute(evalDWRGrad,[2 1 3]));
XZ = matMul(evalGrad,permute(evalBasis,[2 1 3]));
XYXZ = matMul(permute(XY,[2 1 3]),XZ);

if norm(absGrad) > 0
    term1 = matMul(reshape(evalD2W./absGrad.^2,[1 1 length(x)]),XYXZ);
    term2 = -matMul(reshape(evalDW./absGrad.^2,[1 1 length(x)]),XYXZ);
else
    term1 = 0;
    term2 = 0;
end

term3 = matMul(reshape(evalDW,[1 1 length(x)]),YZ); 

D2WTerm = term1 + term2 + term3;
DWTerm = matMul(reshape(evalDW,[1 1 length(x)]),XZ); 

val = 2*DWTerm.*D2WTerm;
