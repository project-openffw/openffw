function p = P3createLinSys(p)

% author: Joscha Gedicke

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Nb = p.level(end).geom.Nb;
Db = p.level(end).geom.Db;

% load problem definition
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
dofU4e = p.level(end).enum.dofU4e;
e4ed = p.level(end).enum.e4ed;
midPoint4e = p.level(end).enum.midPoint4e;
nrElems = p.level(end).nrElems;

% additional data
C = p.statics.basisCoefficients;
curLvl = length(p.level);
degree = loadField('p.params','rhsIntegtrateExactDegree',p,3);

%% Assembling global energy matrix   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = zeros(10,10,nrElems);
Mama = zeros(10,10,nrElems);

genericMama =  ...
 [  420    210    210     84     42     84     14      0    -14     14      
    210    420    210     84     84     42    -14     14      0     14      
    210    210    420     42     84     84      0    -14     14     14      
     84     84     42     28     14     14      0      2     -2      4      
     42     84     84     14     28     14     -2      0      2      4      
     84     42     84     14     14     28      2     -2      0      4      
     14    -14      0      0     -2      2      3     -1     -1      0      
      0     14    -14      2      0     -2     -1      3     -1      0      
    -14      0     14     -2      2      0     -1     -1      3      0      
     14     14     14      4      4      4      0      0      0      1  ];
genericMama = 1/2520 * (C*genericMama*C');

localStima = integrate(n4e,curLvl,degree,@funcHandleStima,p);
localStima = permute(localStima,[2 3 1]); 

localDama = integrate(n4e,curLvl,degree,@funcHandleDama,p);
localDama = permute(localDama,[2 3 1]); 

mu4e = mu(midPoint4e(:,1),midPoint4e(:,2),p);

for curElem = 1:nrElems	

	area = area4e(curElem);
    curMu = mu4e(curElem);
	localMama  = curMu*area*genericMama;
	
	S(:,:,curElem) = localStima(:,:,curElem) + localMama + localDama(:,:,curElem);
    Mama(:,:,curElem) = area*genericMama;
    
end

[I,J] = localDoFtoGlobalDoF(dofU4e,dofU4e);
A = sparse(I(:),J(:),S(:));
B = sparse(I(:),J(:),Mama(:));

%% Assembling Righthandside		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4e = integrate(n4e,curLvl,degree,@funcHandleRHSVolume,p);
b = accumarray(dofU4e(:),f4e(:));

%% Include boundary conditions	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(Nb)  
     g4NbEd = integrate(Nb,curLvl,degree,@funcHandleRHSNb,p);
     n4NbElem = dofU4e(e4ed(NbEd),:);
     neumann = accumarray(n4NbElem(:),g4NbEd(:),[nrNodes+2*nrEdges+nrElems,1]);     
     b = b + neumann;
end
 
u = zeros(nrNodes+2*nrEdges+nrElems,1);
firstPoint4ed = 1/3*( c4n(Db(:,2),:) - c4n(Db(:,1),:) ) + c4n(Db(:,1),:);
secondPoint4ed = 2/3*( c4n(Db(:,2),:) - c4n(Db(:,1),:) )+ c4n(Db(:,1),:);
u(fixedNodes) = u_D([c4n(unique(Db),1);firstPoint4ed(:,1);secondPoint4ed(:,1)],...
              [c4n(unique(Db),2);firstPoint4ed(:,2);secondPoint4ed(:,2)],p);
b = b - A*u;

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A = A;
p.level(end).B = B;
p.level(end).b = b;
p.level(end).x = u;
