function p = TWRTDWRgetFuncVal(p)
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
D2W = p.problem.conjNonLinearFuncSecDer;
Psigma_h = p.statics.sigma_h;     %primal solution X
Dsigma_h = p.statics.DWRsigma_h;  %dual solution   Y
stressBasis = p.statics.stressBasis;     % basis   Z

evalPSigma = Psigma_h(x,y,curElem,lvl,p);
evalDSigma = Dsigma_h(x,y,curElem,lvl,p);
evalBasis = stressBasis(x,y,curElem,lvl,p);

evalD2W = D2W(evalPSigma(:,1),evalPSigma(:,2),curElem,lvl,p);

evalDSigma = reshape(evalDSigma',[1 2 length(x)]);
WQ = matMul(evalDSigma,permute(evalBasis,[2 1 3]));

val = matMul(reshape(evalD2W,[1 1 length(x)]),WQ);

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

%n4e = p.level(end).geom.n4e;
%c4n = p.level(end).geom.c4n;
%curNodes = n4e(curElem,:);
%Coords = c4n(curNodes,:);
stressBasis = p.statics.stressBasis; 
evalBasis = stressBasis(x,y,curElem,lvl,p);
%evalBasis = stressBasis(Coords(:,1),Coords(:,2),curElem,lvl,p);

f4e = p.level(end).DWRf4e;
evalF4e = f4e(curElem);
F4e = evalF4e*ones(2,length(x));
F4e = reshape(F4e,[1 2 length(x)]);

FQ = matMul(F4e,permute(evalBasis,[2 1 3]));

val = FQ;


