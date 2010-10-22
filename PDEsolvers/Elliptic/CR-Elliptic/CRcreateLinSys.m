function p = CRcreateLinSys(p)
%createLinSys.m creates the energy matrix A 
%and the right-hand side b for a nonconforming CR-FE method. 
%The differential operator is full elliptic with piecewise constant
%coefficients (one point gauss integration). 
%
%authors: David Guenther, Jan Reininghaus

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load geometry
Nb = p.level(end).geom.Nb;
n4e = p.level(end).geom.n4e;

% load problem definition
kappa = p.problem.kappa;
mu = p.problem.mu;
lambda = p.problem.lambda;
u_D = p.problem.u_D;
g = p.problem.g;

% load enumerated data
midPoint4ed = p.level(end).enum.midPoint4ed;
midPoint4e = p.level(end).enum.midPoint4e;
DbEd = p.level(end).enum.DbEd;
NbEd = p.level(end).enum.NbEd;
length4ed = p.level(end).enum.length4ed;
gradNC4e = p.level(end).enum.gradNC4e;
ed4e = p.level(end).enum.ed4e;
area4e = p.level(end).enum.area4e;
normals4NbEd = p.level(end).enum.normals4NbEd;
dofU4e = p.level(end).enum.dofU4e;
e4ed = p.level(end).enum.e4ed;

nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;

curLvl = length(p.level);
degree = loadField('p.params','rhsIntegtrateExactDegree',p,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling global energy matrix				   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = zeros(3,3,nrElems);

kappa4e = kappa(midPoint4e(:,1),midPoint4e(:,2),p);
lambda4e = lambda(midPoint4e(:,1),midPoint4e(:,2),p);
mu4e = mu(midPoint4e(:,1),midPoint4e(:,2),p);

genericMama = 1/3 * [1 0 0;
                     0 1 0; 
                     0 0 1];

genericDama = 1/3 * [1 1 1];

for curElem = 1:nrElems
	grad = gradNC4e(:,:,curElem);
	area = area4e(curElem);
	curKappa = kappa4e(:,:,curElem);
    curMu = mu4e(curElem);
    curLambda = lambda4e(curElem,:)';
    
	localStima = grad*curKappa*grad'*area;
	localMama = curMu*genericMama*area;
    localDama = grad*curLambda*area*genericDama;
    
    S(:,:,curElem) = localStima + localDama + localMama;
end

[I,J] = localDoFtoGlobalDoF(dofU4e,dofU4e);
A = sparse(I(:),J(:),S(:));

B = repmat(eye(3)/3,[1,1,nrElems]).*area;
B = sparse(I(:),J(:),B(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling Righthandside					   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f4e = integrate(n4e,curLvl,degree,@funcHandleRHSVolume,p);
% if isempty(f4e)
%     f4e = integrate(n4e,curLvl,degree,@funcHandleRHSVolume,p);
    f4e = integrate(n4e,curLvl,degree,@RHS,p);
% end

b = accumarray(ed4e(:),f4e(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Include boundary conditions						  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if ~isempty(Nb)  
%      g4NbEd = length4ed(NbEd).*...
%              g(midPoint4ed(NbEd,1),midPoint4ed(NbEd,2),normals4NbEd,p);
     g4NbEd = integrate(Nb,curLvl,degree,@funcHandleRHSNb,p);
     ed4NbElem = ed4e(e4ed(NbEd),:);
     neumann = accumarray(ed4NbElem(:),g4NbEd(:),[nrEdges,1]);
     b = b + neumann;
 end

x = zeros(nrEdges,1);
x(DbEd) = u_D(midPoint4ed(DbEd,1),midPoint4ed(DbEd,2),p);
b = b - A*x;

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A = A;
p.level(end).B = B;
p.level(end).b = b;
p.level(end).x = x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = RHS(x,y,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
grad4e = p.level(lvl).enum.gradNC4e;

evalGrad = gradU_exact(x,y,p);
evalBasis = grad4e(:,:,curElem);

val = evalBasis*evalGrad';
val = reshape(val,[3 1 length(x)]);