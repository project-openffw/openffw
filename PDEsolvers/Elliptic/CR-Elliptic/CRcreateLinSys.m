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

f4e = integrate(n4e,curLvl,degree,@funcHandleRHSVolume,p);
% b = accumarray(ed4e(:),f4e(:));
b = full(sparse(ed4e(:),ones(size(ed4e,1)*size(ed4e,2),1),f4e(:),nrEdges,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Include boundary conditions						  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if ~isempty(Nb)  
%      g4NbEd = length4ed(NbEd).*...
%              g(midPoint4ed(NbEd,1),midPoint4ed(NbEd,2),normals4NbEd,p);
     g4NbEd = integrate(Nb,curLvl,degree,@funcHandleRHSNb,p);
     ed4NbElem = ed4e(e4ed(NbEd),:);
%      neumann = accumarray(ed4NbElem(:),g4NbEd(:),[nrEdges,1]);
     neumann = sparse(ed4NbElem(:),ones(size(ed4NbElem,1)*size(ed4NbElem,2),1), ...
                           g4NbEd(:),nrEdges,1);     
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
