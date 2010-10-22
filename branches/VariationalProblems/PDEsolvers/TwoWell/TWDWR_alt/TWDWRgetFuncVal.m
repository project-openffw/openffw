function p = TWDWRgetFuncVal(p)
% author: Lena Noack
%% INPUT
n4e = p.level(end).geom.n4e;

length4ed = p.level(end).enum.length4ed;
ed4e = p.level(end).enum.ed4e;
e4ed = p.level(end).enum.e4ed;
n4ed = p.level(end).enum.n4ed;
DbEd = p.level(end).enum.DbEd;

nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;

x = p.level(end).x0;
lvl = size(p.level,2);

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute the function value E(x)

BT = zeros(1,3,nrElems);
CT = zeros(1,3,nrElems);
DT = zeros(1,1,nrElems);
ET = zeros(1,3,nrElems);
FT = zeros(1,3,nrElems);
GT = zeros(1,1,nrElems);

for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curEdgeLengths = length4ed(curEdges);
    curX4ed = x(curEdges);
    curX4e = x(nrEdges + curElem);
    curLambda1 = x(nrEdges + nrElems + curEdges);
    curLambda2 = x(nrEdges + nrElems + nrEdges + curElem);
    
    signum = ones(1,3);
	I = find(e4ed(curEdges,2) == curElem);
	signum(I) = -1;    
    %divergence of basis functions on T,e.g. div(q) = sign*|E|/|T|
    div_qh = signum.*curEdgeLengths';
    
    %calc int_T  DW(p_h)*q_h
    DW = integrate(n4e(curElem,:),lvl,degree,@integrand1,p);

    %calc int_T  D^2W(p_h)*\lambda1*q_h
    D2W = integrate(n4e(curElem,:),lvl,degree,@integrand2,p);
    D2W = D2W(:)';

    BT(:,:,curElem) = DW;
   	%calc int_T u_h*div(q_h)
    CT(:,:,curElem) = curX4e*div_qh;	
    %calc int_T v_h*div(p_h)
    DT(:,:,curElem) = curX4ed'*div_qh';	

    ET(:,:,curElem) = DW - D2W;
    FT(:,:,curElem) = curLambda2*div_qh;
    
    GT(:,:,curElem) = curLambda1'*div_qh';	
end

I = ed4e';
B = accumarray(I(:),BT(:));
C = accumarray(I(:),CT(:));
E = accumarray(I(:),ET(:));
F = accumarray(I(:),FT(:));

I = 1:nrElems;
D = accumarray(I(:),DT(:));
G = accumarray(I(:),GT(:));

boundary = zeros(nrEdges + nrElems,1);
boundary(DbEd) = integrate(n4ed(DbEd,:),lvl,degree,@integrandBoundary,p);

b = -boundary;
f4e = p.level(end).f4e;
b(nrEdges+1:end) = f4e;

% b = -boundary;
% f4e = p.level(end).f4e;
% b(nrEdges+1:end) = f4e;

funcVal1 = [B + C;
               D]   +   b;

funcVal2 = [E - F;
            G];

funcVal = [funcVal2;
           funcVal1];

%% OUTPUT 
p.level(end).funcVal = funcVal;

%% supply integrand DW(\nabla u_h)*\nabla w_h + 2(u_h - f)w_h
function val = integrand1(x,y,curElem,lvl,p)
u_h = p.statics.u_h;
f = p.problem.f;
evalUh = u_h(x,y,curElem,lvl,p);
evalF = f(x,y,curElem,lvl,p);
basisU = p.statics.basisP1; %which basis?
evalBasisU = basisU(x,y,curElem,lvl,p);
evalBasisU = reshape(evalBasisU,[3 1 length(x)]);
evalDif = reshape(2*(evalUh-evalF)',[1 1 length(x)]);

rhs = matMul(evalDif,evalBasisU);%rhs = 2(u_h-f)w_h


F1 = p.problem.F1;
F2 = p.problem.F2;
CONV = p.params.CONV;

if strcmp(CONV,'c')
W = p.problem.nonLinearReg;
else
W = p.problem.nonLinearExact;
end

sigma_h = p.statics.sigma_h;
evalSigma = sigma_h(x,y,curElem,lvl,p);
evalSigma = reshape(evalSigma',[2 1 length(x)]);
absSigma = ( evalSigma(:,1).^2 + evalSigma(:,2).^2 ).^(1/2);

evalW = W(absSigma(1),absSigma(2),curElem,lvl,p)*evalSigma; % right functions?

size(evalW)
size(rhs)
val = evalW + rhs;

%% supply integrand: D2W(\nabla u_h)*\nabla w_h\nabla q_h + 2 w_h q_h
function val = integrand2(x,y,curElem,lvl,p)
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

sigma_h = p.statics.sigma_h;
evalSigma = sigma_h(x,y,curElem,lvl,p);

absSigma = ( evalSigma(:,1).^2 + evalSigma(:,2).^2 ).^(1/2);
grad1 = evalSigma(:,1);
grad2 = evalSigma(:,2);

stressBasis = p.statics.stressBasis;
basisU = p.statics.basisP1;%which basis?

evalGrad = evalSigma;
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

val = term1 + term2 + 2*wq; size(val)

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