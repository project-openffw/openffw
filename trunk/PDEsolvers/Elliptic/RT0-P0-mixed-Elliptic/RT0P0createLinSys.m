function p = RT0P0createLinSys(p)
% creates the energy matrix A 
% and the right-hand side b for the mixed RT0-P0-FE method. 
% The differential operator is full elliptic with piecewise constant
% coefficients \lambda and \mu (one-point gauss integration). \kappa
% has to be the identity matrix.

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
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;
Nb = p.level(end).geom.Nb;

% load problem definition
g = p.problem.g;
u_D = p.problem.u_D;

mu = p.problem.mu;
% kappa = p.problem.kappa;
lambda = p.problem.lambda;

% load enumerated data
tangents4e = p.level(end).enum.tangents4e;
midPoint4ed = p.level(end).enum.midPoint4ed;
midPoint4e = p.level(end).enum.midPoint4e;
length4ed = p.level(end).enum.length4ed;
e4ed = p.level(end).enum.e4ed;
DbEd = p.level(end).enum.DbEd;
NbEd = p.level(end).enum.NbEd;
ed4e = p.level(end).enum.ed4e;
area4e = p.level(end).enum.area4e;
normals4NbEd = p.level(end).enum.normals4NbEd;
dofU4e = p.level(end).enum.dofU4e;
dofSigma4e = p.level(end).enum.dofSigma4e;

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
BT = zeros(3,3,nrElems);
CT = zeros(1,3,nrElems);
DT = zeros(1,3,nrElems);
FT = zeros(1,3,nrElems);

M = [2 0 1 0 1 0;
 	 0 2 0 1 0 1;
	 1 0 2 0 1 0;
	 0 1 0 2 0 1;
	 1 0 1 0 2 0;
	 0 1 0 1 0 2];

N = zeros(6,3);

kappa4e = repmat([1 0; 0 1],[1,1,nrElems]);
lambda4e = lambda(midPoint4e,p);
mu4e = mu(midPoint4e,p);

for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curEdgeLengths = length4ed(curEdges)';
	curNodes = n4e(curElem,:);
	curTangents = tangents4e(:,:,curElem)' .* [curEdgeLengths;curEdgeLengths];
   	curKappa = kappa4e(:,:,curElem);
	curLambda = lambda4e(curElem,:)';
	curCoords = c4n(curNodes([3 1 2]),:)';
	area = area4e(curElem);
	curMidPoint = repmat(midPoint4e(curElem,:),3,1);
	
	signum = ones(1,3);
	I = find(e4ed(curEdges,2) == curElem);
	signum(I) = -1;
	%divergence of basis functions on T,e.g. div(q) = sign*|E|/|T|
    div_qh = signum.*curEdgeLengths;
    
    
	N(:,1) = [0;0;curTangents(:,3);	-curTangents(:,2)];
	N(:,2) = [-curTangents(:,3);0;0; curTangents(:,1)];
	N(:,3) = [curTangents(:,2);	-curTangents(:,1);0;0];
	L = diag(div_qh);
    	
    %calc int_T  p_h*q_h
    BT(:,:,curElem) = 1/(48*area)* L * N' * M * N * L;
    %calc int_T u_h*div(q_h)
   	CT(:,:,curElem) = div_qh;
    %calc int_T v_h*div(p_h)
    FT(:,:,curElem) = div_qh;
    %calc int_T (\lambda*p_h)*v_h
	DT(:,:,curElem) = 1/2*signum.*curEdgeLengths.*(curLambda'*(curMidPoint'-curCoords));
end

% bilinear form b(p,q)
% local DoF -> global DoF
[I,J] = localDoFtoGlobalDoF(dofSigma4e,dofSigma4e);
B = sparse(I(:),J(:),BT(:));

% bilinear form c(p,v)
[I,J] = localDoFtoGlobalDoF(dofU4e,dofSigma4e);
C = sparse(I(:),J(:),CT(:),nrElems,nrEdges);
F = sparse(I(:),J(:),FT(:),nrElems,nrEdges);

% bilinear form d(p,v)
[I,J] = localDoFtoGlobalDoF(dofU4e,dofSigma4e);
D = sparse(I(:),J(:),DT(:),nrElems,nrEdges);

% bilinear form e(u,v)
%calc int_T u_h*v_h
[I,J] = localDoFtoGlobalDoF(dofU4e,dofU4e);
E = sparse(I,J,mu4e.*area4e);

% global energy matrix
% \int (p_h*q_h + u_h*div(q_h) = \int u_D*(q_h*\nu)
% \int (div(p_h)*v_h - (\lambda*p_h)*v_h - \mu*u_h*v_h) = -\int f*v_h
A = [B			C';
	 (F-D)		-E];

%% Assembling Righthandside %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc \int u_D*(q_h*\nu)
% note that (q_h*\nu)(midPoint4ed) = 1
b = zeros(nrEdges + nrElems,1);
% b(DbEd) = length4ed(DbEd).*u_D(midPoint4ed(DbEd,:),p);

f4e = integrate(n4e,curLvl,degreeRhs,@funcHandleRHSVolume,p);
b(nrEdges+1:end) = -f4e;

%% Include Neumann conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = zeros(nrEdges + nrElems,1);
if ~isempty(Nb)  
	x(NbEd) = g(midPoint4ed(NbEd,:),normals4NbEd,p);
	b = b - A*x;
end

%% Include Dirichlet conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b(DbEd) = length4ed(DbEd).*u_D(midPoint4ed(DbEd,:),p);

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A = A;
p.level(end).b = b;
p.level(end).x = x;
p.level(end).BT = BT;