function p = P1P1createLinSys(p)
%createLinSys.m creates the energy matrix A 
%and the right-hand side b for a P1-FE method 
%in linear elasticity. 
%
%authors: David Guenther, Jan Reininghaus
%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = p.PDE.lambda;
mu = p.PDE.mu;

n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Db = p.level(end).geom.Db;
Nb = p.level(end).geom.Nb;

g = p.problem.g;
u_D = p.problem.u_D;

midPoint4ed = p.level(end).enum.midPoint4ed;
area4e = p.level(end).enum.area4e;
grad4e = p.level(end).enum.grad4e;
normals4NbEd = p.level(end).enum.normals4NbEd;
length4ed = p.level(end).enum.length4ed;
NbEd = p.level(end).enum.NbEd;
dofU4e = p.level(end).enum.dofU4e;

nrElems = p.level(end).nrElems;
nrNodes = p.level(end).nrNodes;

curLvl = length(p.level);
degree = loadField('p.params','rhsIntegtrateExactDegree',p,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = mu*[2,0,0;0,2,0;0,0,1] + lambda*[1,1,0;1,1,0;0,0,0];
R = zeros(3,6);
S = zeros(6,6,nrElems);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling global energy matrix				   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for curElem = 1:nrElems
    area = area4e(curElem);
    grad = grad4e(:,:,curElem);
    
    % P1 approximation of the first component of u 
    R([1,3],[1,2,3]) = grad';
    % P1 approximation of the second component of u
    R([3,2],[4,5,6]) = grad';

    stema = area*R'*C*R;
    S(:,:,curElem) = stema;
end

% local DoF -> global DoF
[I,J] = localDoFtoGlobalDoF(dofU4e,dofU4e);
A = sparse(I(:),J(:),S(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling Righthandside					   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4e = integrate(n4e,curLvl,degree,@funcHandleRHSVolume,p);

intFU = f4e(:,1:3)';
intFV = f4e(:,4:6)';
n4eT = n4e';
I = [n4eT(:); n4eT(:) + nrNodes];
S = [intFU(:); intFV(:)];
b = accumarray(I,S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Include boundary conditions						  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Db)
    DirichletNodes = unique(Db);
    nrDbNodes = length(DirichletNodes);
    
    DbVal = u_D(c4n(DirichletNodes,1),c4n(DirichletNodes,2),p);

    I = 1:nrDbNodes;
    J = DirichletNodes;
    S = ones(nrDbNodes,1);
    BU = sparse(I,J,S,nrDbNodes,2*nrNodes);

    J = nrNodes + DirichletNodes; 
    BV = sparse(I,J,S,nrDbNodes,2*nrNodes);

    B = [BU;BV];
end

if ~isempty(Nb)
    g4Nb = g(midPoint4ed(NbEd,1),midPoint4ed(NbEd,2),normals4NbEd,p);
    intG = g4Nb.*[length4ed(NbEd),length4ed(NbEd)]/2;
    
    bTempU = zeros(size(b));
    bTempV = zeros(size(b));
    
    for k = 1:length(NbEd)
        curNb = Nb(k,:);
        bTempU(curNb) = bTempU(curNb) + intG(k,1)*[1;1];
        bTempV(nrNodes + curNb) = bTempV(nrNodes + curNb) + intG(k,2)*[1;1];
    end

    b = b + bTempU + bTempV;
end

A = [A      B';...
     B      sparse(size(B,1),size(B,1))];
 
b = [b;DbVal(:)];

freeNodes = 1:size(A,1);

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A = A;
p.level(end).b = b;
p.level(end).enum.freeNodes = freeNodes';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
