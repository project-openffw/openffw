function p = P2createLinSys(p)
% creates the energy matrix A 
% and the right-hand side b for a conforming P2-FE method. 
% The differential operator is full elliptic with piecewise constant
% coefficients 

% Copyright 2007 Joscha Gedicke
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
Db = p.level(end).geom.Db;

% load problem definition
kappa = p.problem.kappa;
lambda = p.problem.lambda;
mu = p.problem.mu;
u_D = p.problem.u_D;

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
f4e = loadField('p.level(end)','f4e',p,[]);

% load integration parameters
% degreeStima = p.params.integrationDegrees.createLinSys.Stima;
% degreeDama = p.params.integrationDegrees.createLinSys.Dama;
% degreeMama = p.params.integrationDegrees.createLinSys.Mama;
degreeRhs = p.params.integrationDegrees.createLinSys.Rhs;
degreeNeumann = p.params.integrationDegrees.createLinSys.Neumann;

% get current level number
curLvl = length(p.level);

%% Assembling global energy matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = zeros(6,6,nrElems);

genericMama = 1/360*[  6 -1 -1  0 -4  0 
                      -1  6 -1  0  0 -4
                      -1 -1  6 -4  0  0
                       0  0 -4 32 16 16
                      -4  0  0 16 32 16
                       0 -4  0 16 16 32 ];

kappa4e = kappa(midPoint4e,p);
lambda4e = lambda(midPoint4e,p);
mu4e = mu(midPoint4e,p);

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
end

% localStima = integrateVectorised(n4e,curLvl,degreeStima,@funcHandleStimaVectorised,p);
% localDama = integrateVectorised(n4e,curLvl,degreeDama,@funcHandleDamaVectorised,p);
% localMama = integrateVectorised(n4e,curLvl,degreeMama,@funcHandleMamaVectorised,p);
% 
% localStima = integrate(n4e,curLvl,degreeStima,@funcHandleStima,p);
% localDama = integrate(n4e,curLvl,degreeDama,@funcHandleDama,p);
% localMama = integrate(n4e,curLvl,degreeMama,@funcHandleMama,p);
% 
% S2 = permute(localStima+localDama+localMama,[2 3 1]);

[I,J] = localDoFtoGlobalDoF(dofU4e,dofU4e);
A = sparse(I(:),J(:),S(:));


%% Assembling Righthandside %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~f4e)
    % f4e = integrate(n4e,curLvl,degreeRhs,@funcHandleRHSVolume,p);
    f4e = integrateVectorised(n4e,curLvl,degreeRhs,@funcHandleRHSVolumeVectorised,p);
end

b = accumarray(dofU4e(:),f4e(:));


%% Include Neumann conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Nb)  
%      g4NbEd = integrate(Nb,curLvl,degreeNeumann,@funcHandleRHSNb,p);
     g4NbEd = integrateVectorised(Nb,curLvl,degreeNeumann,@funcHandleRHSNbVectorised,p);
     n4NbElem = dofU4e(e4ed(NbEd),:);
     neumann = accumarray(n4NbElem(:),g4NbEd(:),[nrNodes+nrEdges,1]);     
     b = b + neumann;
end

%% Include Dirichlet conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = zeros(nrNodes+nrEdges,1);
u(fixedNodes) = u_D([c4n(unique(Db),:);midPoint4ed(DbEd,:)],p);
b = b - A*u;

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A = A;
p.level(end).b = b;
p.level(end).x = u;