function p = AWcreateLinSys(p)
%createLinSys.m creates the energy matrix A 
%and the right-hand side b for the Arnold-Winther mixed FE 
%in linear elasticity. 
%
%authors: David Guenther, Jan Reininghaus

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda  =  p.PDE.lambda;
mu  =  p.PDE.mu;

g = p.problem.g;
u_D = p.problem.u_D;

n4e  =  p.level(end).geom.n4e;
c4n  =  p.level(end).geom.c4n;
Db  =  p.level(end).geom.Db;
Nb  =  p.level(end).geom.Nb;

nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;

DbEd = p.level(end).enum.DbEd;
NbEd = p.level(end).enum.NbEd;
e4ed = p.level(end).enum.e4ed;
ed4n = p.level(end).enum.ed4n;
n4ed = p.level(end).enum.n4ed;

area4e = p.level(end).enum.area4e;
angle4n = p.level(end).enum.angle4n;
normals4DbEd = p.level(end).enum.normals4DbEd;
normals4NbEd = p.level(end).enum.normals4NbEd;
tangents4e = p.level(end).enum.tangents4e;
length4ed = p.level(end).enum.length4ed;
dofU4e = p.level(end).enum.dofU4e;
dofSigma4e = p.level(end).enum.dofSigma4e;
grad4e = p.level(end).enum.grad4e;
basisCoefficents = p.level(end).basisCoefficents;
massMatrix = p.level(end).enum.massMatrix;

curLvl = length(p.level);
degree = loadField('p.params','rhsIntegtrateExactDegree',p,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dof4e = [dofSigma4e,dofU4e];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling global energy matrix				   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%entries of the tensor Cinv
lame1 = (lambda + 2*mu)/(4*mu*(lambda+mu));
lame2 = -lambda /(4*mu*(lambda+mu));
lame3 = 1/mu;

I = (1:3:28);

% coeff a^(j)
coeffA = basisCoefficents(I,:,:);
% coeff b^(j)
coeffB = basisCoefficents(I+1,:,:);
% coeff c^(j)
coeffC = basisCoefficents(I+2,:,:);

genericMatrix1 = [ ones(3,3),       zeros(3,3);
    zeros(3,3),       ones(3,3)   ] + eye(6,6);

genericMatrix2 = eye(6);
genericMatrix2 = genericMatrix2([1,3,5,2,4,6],:);

L = zeros(2,3,3);
dummyZeros = zeros(2,3);

S = zeros(30,30,nrElems);

matrixD = zeros(6,24,nrElems);

for curElem = 1:nrElems
    area = area4e(curElem);
    curGrad = grad4e(:,:,curElem);
    curBasisCoefficents = basisCoefficents(:,:,curElem);

    %%%%%%%%Construction of matrix A^1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    curCoeffA = coeffA(:,:,curElem);
    curCoeffB = coeffB(:,:,curElem);
    curCoeffC = coeffC(:,:,curElem);

    A1 = lame1 * curCoeffA' * massMatrix * curCoeffA;
    A2 = lame1 * curCoeffB' * massMatrix * curCoeffB;
    A3 = lame2 * curCoeffA' * massMatrix * curCoeffB;
    A4 = lame2 * curCoeffB' * massMatrix * curCoeffA;
    A5 = lame3 * curCoeffC' * massMatrix * curCoeffC;

    A = area * (A1 + A2 + A3 + A4 + A5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%Construction of matrix A^2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:3
        L(:,:,j) = [curGrad(j,1), 0             , curGrad(j,2);
            0           , curGrad(j,2)  , curGrad(j,1)];
    end

    L1 = repmat([L(:,:,1),L(:,:,2),L(:,:,3)],3,1);

    L2 = [	L(:,:,2),	dummyZeros,     L(:,:,3);
        L(:,:,1),	 L(:,:,3),   dummyZeros;
        dummyZeros,	 L(:,:,2),     L(:,:,1)];

    L3 = [  L(:,:,2), dummyZeros,  -L(:,:,3);
        -L(:,:,1),  L(:,:,3),   dummyZeros;
        dummyZeros, -L(:,:,2),   L(:,:,1)];


    L4 = zeros(6,3);

    % dim(D): 6x30 * 30x24
    D = [L1, L2 ,L3, L4] * curBasisCoefficents;

    matrixD(:,:,curElem) = D;

    B = area/12 *genericMatrix1 * genericMatrix2* D;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    B = B*sqrt(area);
    S(:,:,curElem) = [  A     B';
                        B    zeros(6,6)];
end

%local DoF to global DoF
[I,J] = localDoFtoGlobalDoF(dof4e,dof4e);
A = sparse(I(:),J(:),S(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assembling Righthandside                   	      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b = zeros(size(A,1),1);

f4e = integrate(n4e,curLvl,degree,@funcHandleRHSVolume,p);

for curElem  =  1:nrElems
     area = area4e(curElem);
    fs = f4e(curElem,:);
    
    I4 = 3*nrNodes + 4*nrEdges + 3*nrElems + 6*(curElem-1) + (1:6);
    b(I4) = -fs * sqrt(area);
%     b(I4) = -fs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Include boundary conditions						  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Dirichlet boundary
for k  =  1:length(DbEd)
    curDbEd = DbEd(k);
    curEdgeLength = length4ed(curDbEd);
    curElem = e4ed(curDbEd,1);
    curNodes = Db(k,:);
    curCoords = c4n(curNodes(1),:);
    curTangents = tangents4e(:,:,curElem);
    curBasisCoefficents = basisCoefficents(:,:,curElem);
    n = normals4DbEd(k,:);
    rhs = zeros(30,1);

    n1 = find(Db(k,1) == n4e(curElem,:));
    n2 = rem(n1,3) + 1;

    val = quadv(@intud,0,1,[],[],curCoords,curEdgeLength,...
        curTangents(n1,:),u_D,p);

    In1_1 = [20 0];
    In1_2 = [3 6 20 0 23 26];
    In2_1 = rem(In1_1+20,30);
    In2_2 = rem(In1_2+20,30);

    for i = 1:24
        a = curBasisCoefficents(:,i)';
        a = a([1:3:30,2:3:30,3:3:30]);
        rhs(i) = (n*[a(In1_1+n2),a(In1_2+n1);
            a(In2_1+n2),a(In2_2+n1)]) * val;
    end

    globalDoFs = dof4e(curElem,:);
    b(globalDoFs) = b(globalDoFs) + rhs;
end

%Neumann boundary
if ~isempty(NbEd)
    nrNeumannDoFs = 4;
    b4Nb = zeros(nrNeumannDoFs*length(NbEd),1);
    dof4Nb = zeros(nrNeumannDoFs,length(NbEd));
    
    %include neumann data for moment zero and moment one
    for k = 1:length(NbEd)
        curNbEd = NbEd(k);
        curEdgeLength = length4ed(curNbEd);
        curElem = e4ed(curNbEd,1);
        curNodes = Nb(k,:);
        curCoords = c4n(curNodes(1),:);
        curTangents = tangents4e(:,:,curElem);
        curNormal = normals4NbEd(k,:);

        n1 = find(Nb(k,1) == n4e(curElem,:));

        localIndex4MomentZero = [10 11] + 4*(n1-1);
        localIndex4MomentOne  = [12 13] + 4*(n1-1);
        localIndex = [localIndex4MomentZero,localIndex4MomentOne];
        dof4Nb(:,k) = dof4e(curElem,localIndex);

        b4Nb(4*(k-1)+1:4*k) = curEdgeLength^2*quadv(@intg,0,1,1e-10,[],curCoords,curEdgeLength,...
                                    curTangents(n1,:),curNormal,g,p);
                                
%                 b4Nb(4*(k-1)+1:4*k) = quadv(@intg,0,1,1e-10,[],curCoords,curEdgeLength,...
%                                     curTangents(n1,:),curNormal,g,p);
    end

    I4Nb = 1:nrNeumannDoFs*length(NbEd);
    NbEdLengths = length4ed(NbEd).^2;
    S =  repmat(NbEdLengths',4,1);  
%     S = ones(nrNeumannDoFs,length(NbEd));
    
    V = sparse(I4Nb(:),dof4Nb(:),S(:),nrNeumannDoFs*length(NbEd),size(A,1));

    A = [   A, V';
            V, sparse(nrNeumannDoFs*length(NbEd),nrNeumannDoFs*length(NbEd))];

    b = [   b;
        b4Nb];
    
    %include neumann data for the neumann nodes
    NbNodes = unique(Nb);
    normals4NbNodes = zeros(2,2,nrNodes);

    for k = 1:length(NbEd)
        curNbEdge = NbEd(k);
        curNormal = normals4NbEd(k,:);
        curNodes = n4ed(curNbEdge,:);
        normals4NbNodes(2,:,curNodes) =  normals4NbNodes(1,:,curNodes);
        normals4NbNodes(1,:,curNodes) = repmat(curNormal,2,1)';
    end

    I = Nb(:);
    S = ones(size(I));
    I = accumarray(I,S);
    [S,I] = find(I == 1);

    NbAngles = angle4n(NbNodes);
    cornerNodes = find(abs(NbAngles - pi) > 100*eps);
    cornerNodes = NbNodes(cornerNodes);

    cornerNodes = setdiff(cornerNodes,S);
    straightNodes =  setdiff(NbNodes,cornerNodes);

    if ~isempty(cornerNodes)
        I = zeros(3,3,length(cornerNodes));
        J = zeros(3,3,length(cornerNodes));
        S = zeros(3,3,length(cornerNodes));
        B = zeros(3,length(cornerNodes));
        genericMatrix = repmat([1,2,3],3,1);
        
        for k = 1:length(cornerNodes)
            [curElem,localNode] = find(n4e == cornerNodes(k));
            curDoF = dofSigma4e(curElem(1),:);
            dof = curDoF(3*(localNode(1)-1) + (1:3))';
            I(:,:,k) = genericMatrix + 3*(k-1);
            J(:,:,k) = [dof,dof,dof];
            n = normals4NbNodes(:,:,cornerNodes(k));
            V = [n(1,1) 0 n(1,2);
                0       n(1,2) n(1,1);
                n(2,1) 0 n(2,2);
                0       n(2,2) n(2,1)];
            S(:,:,k) = (V' * V)';
            
            
            curNode = cornerNodes(k);
            % this is just one arbitrary edge
            curNbEd = find(ed4n(curNode,:));
            curNbEd = ed4n(curNode,curNbEd(1));
            curEdgeLength = length4ed(curNbEd);
            S(:,:,k) = S(:,:,k) * curEdgeLength^2;
            
            curCoords = c4n(curNode,:);
            g1 = g(curCoords(1),curCoords(2),n(1,:),p)';
            g2 = g(curCoords(1),curCoords(2),n(2,:),p)';
            
%              G = [g1;g2];
            G = curEdgeLength^2*[g1;g2];
            B(:,k) = V' * G;
        end
        V = sparse(I(:),J(:),S(:),3*length(cornerNodes),size(A,1));

        A = [A, V';
            V, sparse(3*length(cornerNodes),3*length(cornerNodes))];

        b = [b;B(:)];
    end

    if ~isempty(straightNodes)
        I = zeros(3,2,length(straightNodes));
        J = zeros(3,2,length(straightNodes));
        S = zeros(3,2,length(straightNodes));
        B = zeros(2,length(straightNodes));
        genericMatrix = repmat([1,2],3,1);
        
        for k = 1:length(straightNodes)
            [curElem,localNode] = find(n4e == straightNodes(k));
            curDoF = dofSigma4e(curElem(1),:);
            dof = curDoF(3*(localNode(1)-1) + (1:3))';
            I(:,:,k) = genericMatrix + 2*(k-1);
            J(:,:,k) = [dof,dof];
            n = normals4NbNodes(:,:,straightNodes(k));
            V = [n(1,1) 0 n(1,2);
                0       n(1,2) n(1,1)];
            S(:,:,k) = V';
            curNode = straightNodes(k);
            
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5   
            % this is just one arbitrary edge
            curNbEd = find(ed4n(curNode,:));
            curNbEd = ed4n(curNode,curNbEd(1));
            curEdgeLength = length4ed(curNbEd);
            S(:,:,k) = S(:,:,k) * curEdgeLength^2;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            
            curCoords = c4n(curNode,:);
            G = curEdgeLength^2*g(curCoords(1),curCoords(2),n(1,:),p)';
%               G = g(curCoords(1),curCoords(2),n(1,:),p)';
            B(:,k) = G;
        end
        V = sparse(I(:),J(:),S(:),2*length(straightNodes),size(A,1));

        A = [A, V';
            V, sparse(2*length(straightNodes),2*length(straightNodes))];

        b = [b;B(:)];
    end
end

freeNodes = 1:size(A,1);

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A  =  A;
p.level(end).matrixD  =  matrixD;
p.level(end).b  =  b;
p.level(end).enum.freeNodes  =  freeNodes';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = intud(x,x0,edgeLength,tangent,func,p)

%integration along the Dirichlet edge
E = edgeLength*tangent;
x1 = x0(1) + E(1)*x;
y1 = x0(2) + E(2)*x;

DbVal = func(x1,y1,p);

val = [	DbVal(:,2)'.*x;
    DbVal(:,1)'.*x;
    DbVal(:,1)'.*x.*(1-x);
    DbVal(:,1)'.*x.*(1-x).*(1-2*x);
    DbVal(:,2)'.*(1-x);
    DbVal(:,1)'.*(1-x);
    DbVal(:,2)'.*x.*(1-x);
    DbVal(:,2)'.*x.*(1-x).*(1-2*x)];

val = edgeLength*val;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = intg(x,x0,edgeLength,tangent,normal,func,p)

%integration along the Neumann edge
E = edgeLength*tangent;
x1 = x0(1) + E(1)*x;
y1 = x0(2) + E(2)*x;

w = func(x1,y1,normal,p)';
momentZero  = 1/edgeLength;
momentOne   = (x - 0.5)/edgeLength;

val = [ w(1,:).*momentZero;
    w(2,:).*momentZero;
    w(1,:).*momentOne;
    w(2,:).*momentOne];

val = edgeLength*val;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
