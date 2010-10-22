function p = P2createLinSys(p)

% author: Joscha Gedicke

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Nb = p.level(end).geom.Nb;
Db = p.level(end).geom.Db;

% load problem definition
kappa = p.problem.kappa;
lambda = p.problem.lambda;
mu = p.problem.mu;
u_D = p.problem.u_D;
g = p.problem.g;

% load enumerated data
fixedNodes = p.level(end).enum.fixedNodes;
P1grad4e = p.level(end).enum.P1grad4e;
area4e = p.level(end).enum.area4e;
nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;
NbEd = p.level(end).enum.NbEd;
DbEd = p.level(end).enum.DbEd;
dofU4e = p.level(end).enum.dofU4e;
e4ed = p.level(end).enum.e4ed;
midPoint4e = p.level(end).enum.midPoint4e;
midPoint4ed = p.level(end).enum.midPoint4ed;

% additional data
C = p.statics.basisCoefficients;
curLvl = length(p.level);
degree = loadField('p.params','rhsIntegtrateExactDegree',p,1);



%% Assembling global energy matrix	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = zeros(6,6,nrElems);
Mama = zeros(6,6,nrElems);

genericMama = 1/360*[  6 -1 -1  0 -4  0 
                      -1  6 -1  0  0 -4
                      -1 -1  6 -4  0  0
                       0  0 -4 32 16 16
                      -4  0  0 16 32 16
                       0 -4  0 16 16 32 ];

kappa4e = kappa(midPoint4e(:,1),midPoint4e(:,2),p);
lambda4e = lambda(midPoint4e(:,1),midPoint4e(:,2),p);
mu4e = mu(midPoint4e(:,1),midPoint4e(:,2),p);

for curElem = 1:nrElems	
    
	P1grad = P1grad4e(:,:,curElem);
	area = area4e(curElem);
	curKappa = kappa4e(:,:,curElem);
    curLambda = lambda4e(curElem,:)';
    curMu = mu4e(curElem);
    
    A11 = area*P1grad*curKappa*P1grad';
    A12 = 1/3*area*P1grad*curKappa*(P1grad([ 1 2 1],:)+P1grad([ 2 3 3],:))';
    A21 = 1/3*area*(P1grad([ 1 2 1],:)+P1grad([ 2 3 3],:))*curKappa*P1grad';
    A22 = area*1/12*[2,1,1; 1,2,2; 1,2,2].*((P1grad([ 1 2 1],:))*curKappa*(P1grad([ 1 2 1],:))')...
          +area*1/12*[1 2 1; 1 1 1; 1 1 1].*((P1grad([ 1 2 1],:))*curKappa*(P1grad([ 2 3 3],:))')...
           +area*1/12*[1 1 1; 2 1 1; 1 1 1].*((P1grad([ 2 3 3],:))*curKappa*(P1grad([ 1 2 1],:))')...
           +area*1/12*[2 1 2; 1 2 1; 2 1 2].*((P1grad([ 2 3 3],:))*curKappa*(P1grad([ 2 3 3],:))');
	localStima = C*[A11,A12; A21,A22]*C';
    
    B11 = P1grad*curLambda*area*[1 1 1]/3;
    B12 = P1grad*curLambda*area*[1 1 1]/12;
    B21 = P1grad([ 1 2 1],:)*curLambda*area*[1 1 2]/12 + ...
           P1grad([ 2 3 3],:)*curLambda*area*[ 2 2 1]/12;
    B22 = P1grad([ 1 2 1],:)*curLambda*area*[1 1 1]/60 + ...
           P1grad([ 2 3 3],:)*curLambda*area*[ 1 1 1]/60;   
	localDama  = C*[B11,B12; B21,B22]*C';
    
	localMama  = curMu*2*area*genericMama;

    S(:,:,curElem) = localStima + localMama + localDama;
    Mama(:,:,curElem) = 2*area*genericMama;
end

[I,J] = localDoFtoGlobalDoF(dofU4e,dofU4e);
A = sparse(I(:),J(:),S(:));

B = sparse(I(:),J(:),Mama(:));


%% Assembling Righthandside    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f4e = loadField('p.level(end)','f4e',p,[]);
if isempty(f4e)
    f4e = integrate(n4e,curLvl,degree,@funcHandleRHSVolume,p);
end
b = accumarray(dofU4e(:),f4e(:));


%% Include boundary conditions		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(Nb)  
     g4NbEd = integrate(Nb,curLvl,degree+1,@funcHandleRHSNb,p);
     n4NbElem = dofU4e(e4ed(NbEd),:);
     neumann = accumarray(n4NbElem(:),g4NbEd(:),[nrNodes+nrEdges,1]);     
     b = b + neumann;
end
 
u = zeros(nrNodes+nrEdges,1);
u(fixedNodes) = u_D([c4n(unique(Db),1);midPoint4ed(DbEd,1)],...
                    [c4n(unique(Db),2);midPoint4ed(DbEd,2)],p);
b = b - A*u;

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A = A;
p.level(end).B = B;
p.level(end).b = b;
p.level(end).x = u;

