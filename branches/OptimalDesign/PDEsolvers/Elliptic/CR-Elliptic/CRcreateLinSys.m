function p = CRcreateLinSys(p)
% creates the energy matrix A 
% and the right-hand side b for a nonconforming CR-FE method. 
% The differential operator is full elliptic with piecewise constant
% coefficients (one point gauss integration). 

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
% load geometry
Nb = p.level(end).geom.Nb;
n4e = p.level(end).geom.n4e;

% load problem definition
kappa = p.problem.kappa;
mu = p.problem.mu;
lambda = p.problem.lambda;
u_D = p.problem.u_D;

% load enumerated data
midPoint4ed = p.level(end).enum.midPoint4ed;
midPoint4e = p.level(end).enum.midPoint4e;
DbEd = p.level(end).enum.DbEd;
NbEd = p.level(end).enum.NbEd;
gradNC4e = p.level(end).enum.gradNC4e;
ed4e = p.level(end).enum.ed4e;
area4e = p.level(end).enum.area4e;
dofU4e = p.level(end).enum.dofU4e;
e4ed = p.level(end).enum.e4ed;

nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;

% load integration parameters
degreeStima = p.params.integrationDegrees.createLinSys.Stima;
degreeDama = p.params.integrationDegrees.createLinSys.Dama;
degreeMama = p.params.integrationDegrees.createLinSys.Mama;
degreeRhs = p.params.integrationDegrees.createLinSys.Rhs;
degreeNeumann = p.params.integrationDegrees.createLinSys.Neumann;

% get current level number
curLvl = length(p.level);

%% Assembling global energy matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = zeros(3,3,nrElems);

kappa4e = kappa(midPoint4e,p);
lambda4e = lambda(midPoint4e,p);
mu4e = mu(midPoint4e,p);

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

%% Assembling Righthandside	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4e = integrate(n4e,curLvl,degreeRhs,@funcHandleRHSVolume,p);
b = accumarray(ed4e(:),f4e(:));

%% Include Neumann conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Nb)
    g4NbEd = integrate(Nb,curLvl,degreeNeumann,@funcHandleRHSNb,p);
    ed4NbElem = ed4e(e4ed(NbEd),:);
    neumann = accumarray(ed4NbElem(:),g4NbEd(:),[nrEdges,1]);
    b = b + neumann;
end

%% Include Dirichlet conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros(nrEdges,1);
x(DbEd) = u_D(midPoint4ed(DbEd,:),p);
b = b - A*x;

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A = A;
p.level(end).B = B;
p.level(end).b = b;
p.level(end).x = x;