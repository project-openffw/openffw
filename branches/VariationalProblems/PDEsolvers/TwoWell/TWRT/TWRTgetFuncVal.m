function p = TWRTgetFuncVal(p)
% author: David Guenther, Lena Noack 

%% INPUT
n4e = p.level(end).geom.n4e;
length4ed = p.level(end).enum.length4ed;
ed4e = p.level(end).enum.ed4e;
e4ed = p.level(end).enum.e4ed;
n4ed = p.level(end).enum.n4ed;
DbEd = p.level(end).enum.DbEd;
dofU4e = p.level(end).enum.dofU4e;
dofSigma4e = p.level(end).enum.dofSigma4e;

nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;

x = p.level(end).x0;
lvl = size(p.level,2);

f04e = p.level(end).f04e;

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute the function value E(x)

BT = zeros(1,3,nrElems);
CT = zeros(1,3,nrElems);
DT = zeros(1,1,nrElems);
ET = zeros(1,1,nrElems);

for curElem = 1:nrElems    
	curEdges = ed4e(curElem,:);
	curEdgeLengths = length4ed(curEdges);
    curX4ed = x(curEdges);
    curX4e = x(nrEdges + curElem);
    
    signum = ones(1,3);
	signum(e4ed(curEdges,2) == curElem) = -1;    
    %divergence of basis functions on T,e.g. div(q) = sign*|E|/|T|
    div_qh = signum.*curEdgeLengths';
    
    %calc int_T  DW*(p_h)*q_h
    DW = integrate(n4e(curElem,:),lvl,degree,@integrand,p);
    
    BT(:,:,curElem) = DW;
   	%calc int_T u_h*div(q_h)
    CT(:,:,curElem) = curX4e*div_qh;	
    %calc int_T v_h*div(p_h)
    DT(:,:,curElem) = curX4ed'*div_qh';	
    %calc  -2 int_T (u_h-f0)*v_h
    ET(:,:,curElem) = -2*(curX4e-f04e(curElem))*curX4e';	
end

I = dofSigma4e';
B = accumarray(I(:),BT(:));
C = accumarray(I(:),CT(:));
D = accumarray(dofU4e(:),DT(:));
E = accumarray(dofU4e(:),ET(:));

boundary = zeros(nrEdges + nrElems,1);
boundary(DbEd) = integrate(n4ed(DbEd,:),lvl,degree,@integrandBoundary,p);

b = -boundary;
f4e = p.level(end).f4e;
b(nrEdges+1:end) = f4e;

funcVal = [B + C;
               D + E]   +   b;

%% OUTPUT 
p.level(end).funcVal = funcVal;

%% supply integrand DW^*(sigma_h)*q_h = Du_h*q_h
function val = integrand(x,y,curElem,lvl,p)
% W'(|X|)/|X|*X*Y

grad_h = p.statics.grad_h;
stressBasis = p.statics.stressBasis;

evalGrad = grad_h(x,y,curElem,lvl,p);
evalGrad = reshape(evalGrad',[2 1 length(x)]);

evalBasis = stressBasis(x,y,curElem,lvl,p);

val = matMul(evalBasis,evalGrad);

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