function p = P1CRcreateLinSys(p)
% creates the energy matrix A 
% and the right-hand side b for the Kouhia-Stenberg FE 
% in linear elasticity.

% Copyright 2007 Jan Reininghaus, David Guenther
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load problem definition
g = p.problem.g;
u_D = p.problem.u_D;
lambda = p.PDE.lambda;
mu = p.PDE.mu;

% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Db = p.level(end).geom.Db;
Nb = p.level(end).geom.Nb;

% load enumerated data
ed4e = p.level(end).enum.ed4e;
midPoint4ed = p.level(end).enum.midPoint4ed;
normals4NbEd = p.level(end).enum.normals4NbEd;
length4ed = p.level(end).enum.length4ed;
DbEd = p.level(end).enum.DbEd;
NbEd = p.level(end).enum.NbEd;
area4e = p.level(end).enum.area4e;
grad4e = p.level(end).enum.grad4e;
gradNC4e = p.level(end).enum.gradNC4e;
dofU4e = p.level(end).enum.dofU4e;
nrElems = p.level(end).nrElems;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;

% load integration parameters
degreeRhs = p.params.integrationDegrees.createLinSys.Rhs;

% get current level number
curLvl = length(p.level);

% additional data
C = mu*[2,0,0;0,2,0;0,0,1] + lambda*[1,1,0;1,1,0;0,0,0];
R = zeros(3,6);
S = zeros(6,6,nrElems);

%% Assembling global energy matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for curElem = 1:nrElems
    area = area4e(curElem);
    grad = grad4e(:,:,curElem);
    gradNC = gradNC4e(:,:,curElem);

    % P1 approximation of the first component of u 
    R([1,3],[1,2,3]) = grad';
    % CR approximation of the second component of u
    R([3,2],[4,5,6]) = gradNC';

    stema = area*R'*C*R;
    S(:,:,curElem) = stema;
end

% local DoF -> global DoF
[I,J] = localDoFtoGlobalDoF(dofU4e,dofU4e);
A = sparse(I(:),J(:),S(:));


%% Assembling Righthandside %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4e = integrate(n4e,curLvl,degreeRhs,@funcHandleRHSVolume,p);

intFU = f4e(:,1:3)';
intFV = f4e(:,4:6)';
n4eT = n4e';
ed4eT = ed4e';
I = [n4eT(:); ed4eT(:) + nrNodes];
S = [intFU(:); intFV(:)];
b = accumarray(I,S);

%% Include Neumann conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(Nb)
    g4Nb = g(midPoint4ed(NbEd,:),normals4NbEd,p);
    intG = g4Nb.*[length4ed(NbEd),length4ed(NbEd)]/2;

    bTemp = zeros(size(b));
    for k = 1:length(NbEd)
        curNb = Nb(k,:);
        bTemp(curNb) = bTemp(curNb) + intG(k,1)*[1;1];
    end
    b = b + bTemp;

    intGV = 2*intG(:,2);
    I = nrNodes + NbEd;
    b(I) = b(I) + intGV;
end

%% Include Neumann conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(Db)
    DirichletNodes = unique(Db);
    midPoint4Db = midPoint4ed(DbEd,:);
    
    nrDbNodes = length(DirichletNodes);
    nrDbEdges = length(DbEd);

    DbValNodes = u_D(c4n(DirichletNodes,:),p);
    DbValEdges = u_D(midPoint4Db,p);

    B = sparse( nrDbNodes + nrDbEdges,nrNodes + nrEdges );

    for k = 1:nrDbNodes
        curDn = DirichletNodes(k);
        B(k,curDn) = 1;
    end

    for k = 1:nrDbEdges
        curDEd = nrNodes + DbEd(k);
        B(nrDbNodes+k,curDEd) = 1;
    end

    rhs = zeros(nrDbNodes + nrDbEdges,1);
    rhs(1:nrDbNodes) = DbValNodes(:,1);
    rhs(nrDbNodes+1 : end) = DbValEdges(:,2);

    A = [A B';...
         B sparse(size(B,1),size(B,1))];
    b = [b;rhs];
end

x = zeros(size(A,1),1);
freeNodes = 1:size(A,1);

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A = A;
p.level(end).b = b;
p.level(end).x = x;
p.level(end).enum.freeNodes = freeNodes';