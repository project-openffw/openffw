function p = P1createLinSys(p)
%createLinSys.m creates the energy matrix A 
%and the right-hand side b for a conforming P1-FE method. 
%The differential operator is full elliptic with piecewise constant
%coefficients (one point gauss integration). 
%
%authors: David Guenther, Jan Reininghaus

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Nb = p.level(end).geom.Nb;

% load problem definition
kappa = p.problem.kappa;
lambda = p.problem.lambda;
mu = p.problem.mu;
u_D = p.problem.u_D;
g = p.problem.g;

% load enumerated data
fixedNodes = p.level(end).enum.fixedNodes;
midPoint4e = p.level(end).enum.midPoint4e;
midPoint4ed = p.level(end).enum.midPoint4ed;
grad4e = p.level(end).enum.grad4e;
area4e = p.level(end).enum.area4e;
nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
NbEd = p.level(end).enum.NbEd;
length4ed = p.level(end).enum.length4ed;
normals4NbEd = p.level(end).enum.normals4NbEd;
dofU4e = p.level(end).enum.dofU4e;
e4ed = p.level(end).enum.e4ed;
NbEd = p.level(end).enum.NbEd;

curLvl = length(p.level);
degree = loadField('p.params','rhsIntegtrateExactDegree',p,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling global energy matrix				   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = zeros(3,3,nrElems);
genericMama = (ones(3)+eye(3))/12;
genericDama = [1 1 1]/3;

kappa4e = kappa(midPoint4e(:,1),midPoint4e(:,2),p);
lambda4e = lambda(midPoint4e(:,1),midPoint4e(:,2),p);
mu4e = mu(midPoint4e(:,1),midPoint4e(:,2),p);

for curElem = 1:nrElems	
	grad = grad4e(:,:,curElem);
	area = area4e(curElem);
	curKappa = kappa4e(:,:,curElem);
    curLambda = lambda4e(curElem,:)';
    curMu = mu4e(curElem);
    
	localStima = grad*curKappa*grad'*area;
	localDama  = grad*curLambda*area*genericDama;
	localMama  = curMu*area*genericMama;
	
	S(:,:,curElem) = localStima + localMama + localDama;
end

% area = reshape(repmat(area4e',[9 1]),[3 3 nrElems]);
% lambda4e = reshape(lambda4e',[2,1,nrElems]);
% mu4e = reshape(repmat(mu4e',[9 1]),[3 3 nrElems]);
% 
% grad4eT = permute(grad4e,[2 1 3]);
% dama = ones(1,3,nrElems)/3;
mama = repmat((ones(3)+eye(3))/12,[1,1,nrElems]);
% 
% localStima = matMul(matMul(grad4e,kappa4e),grad4eT);
% localDama = matMul(matMul(grad4e,lambda4e),dama);
% localMama = mu4e.*mama;
% 
% S = (localStima + localMama + localDama).*area;

[I,J] = localDoFtoGlobalDoF(dofU4e,dofU4e);
A = sparse(I(:),J(:),S(:));

S = mama.*area;
B = sparse(I(:),J(:),S(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling Righthandside					   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f4e = loadField('p.level(end)','f4e',p,[]);
if isempty(f4e)
    f4e = integrate(n4e,curLvl,degree,@funcHandleRHSVolume,p);
%     f4e = integrate(n4e,curLvl,degree,@RHS,p);
    b = accumarray(n4e(:),f4e(:));
else
    b = accumarray(n4e(:),f4e(:));
%     f4ed = p.level(end).f4ed;
%     I = find(e4ed(:,2)==0);
%     e4ed(I,2) = e4ed(I,1);
%     f4ed =  [f4ed;f4ed];
%     dof = dofU4e(e4ed(:),:);
%     edgePart = accumarray(dof(:),f4ed(:),[nrNodes,1]);     
%     b = b + edgePart;
end

% b = accumarray(n4e(:),f4e(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Include boundary conditions						  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Nb)  
%      g4NbEdalt = length4ed(NbEd)/2.*...
%              g(midPoint4ed(NbEd,1),midPoint4ed(NbEd,2),normals4NbEd,p);
     g4NbEd = integrate(Nb,curLvl,degree,@funcHandleRHSNb,p);
     n4NbElem = n4e(e4ed(NbEd),:);
     neumann = accumarray(n4NbElem(:),g4NbEd(:),[nrNodes,1]);     
     b = b + neumann;
end
 
u = zeros(nrNodes,1);
u(fixedNodes) = u_D(c4n(fixedNodes,1),c4n(fixedNodes,2),p);
b = b - A*u;

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A = A;
p.level(end).B = B;
p.level(end).b = b;
p.level(end).x = u;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = RHS(x,y,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
grad4e = p.level(lvl).enum.grad4e;

evalGrad = gradU_exact(x,y,p);
evalBasis = grad4e(:,:,curElem);

val = evalBasis*evalGrad';
val = reshape(val,[3 1 length(x)]);