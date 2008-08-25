function p = P3createLinSys(p)
% creates the energy matrix A 
% and the right-hand side b for a conforming P3-FE method. 
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
mu = p.problem.mu;
lambda = p.problem.lambda;
u_D = p.problem.u_D;

% load enumerated data
fixedNodes = p.level(end).enum.fixedNodes;
area4e = p.level(end).enum.area4e;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
NbEd = p.level(end).enum.NbEd;
dofU4e = p.level(end).enum.dofU4e;
e4ed = p.level(end).enum.e4ed;
midPoint4e = p.level(end).enum.midPoint4e;
nrElems = p.level(end).nrElems;

% additional data
% C = p.statics.basisCoefficients;

% load integration parameters
degreeStima = p.params.integrationDegrees.createLinSys.Stima;
degreeDama = p.params.integrationDegrees.createLinSys.Dama;
degreeMama = p.params.integrationDegrees.createLinSys.Mama;
degreeRhs = p.params.integrationDegrees.createLinSys.Rhs;
degreeNeumann = p.params.integrationDegrees.createLinSys.Neumann;

% get current level number
curLvl = length(p.level);

%% Assembling global energy matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% genericMama =  ...
%  [  420    210    210     84     42     84     14      0    -14     14      
%     210    420    210     84     84     42    -14     14      0     14      
%     210    210    420     42     84     84      0    -14     14     14      
%      84     84     42     28     14     14      0      2     -2      4      
%      42     84     84     14     28     14     -2      0      2      4      
%      84     42     84     14     14     28      2     -2      0      4      
%      14    -14      0      0     -2      2      3     -1     -1      0      
%       0     14    -14      2      0     -2     -1      3     -1      0      
%     -14      0     14     -2      2      0     -1     -1      3      0      
%      14     14     14      4      4      4      0      0      0      1  ];
% genericMama = 1/2520 * (C*genericMama*C');

localStima = integrateVectorised(n4e,curLvl,degreeStima,@funcHandleStimaVectorised,p);
S = permute(localStima,[2 3 1]); 

lambda4e = lambda(midPoint4e,p);
if nnz(lambda4e) ~= 0
    localDama = integrateVectorised(n4e,curLvl,degreeDama,@funcHandleDamaVectorised,p);
    S = S + permute(localDama,[2 3 1]);
end

mu4e = mu(midPoint4e,p);
if nnz(mu4e) ~= 0
    localMama = integrateVectorised(n4e,curLvl,degreeMama,@funcHandleMamaVectorised,p);
    S = S + permute(localMama,[2 3 1]);
end

[I,J] = localDoFtoGlobalDoF(dofU4e,dofU4e);
A = sparse(I(:),J(:),S(:));


%% Assembling Righthandside %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f4e = integrate(n4e,curLvl,degreeRhs,@funcHandleRHSVolume,p);
f4e = integrateVectorised(n4e,curLvl,degreeRhs,@funcHandleRHSVolumeVectorised,p);
b = accumarray(dofU4e(:),f4e(:));

%% Include Neumann conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(Nb)  
%      g4NbEd = integrate(Nb,curLvl,degreeNeumann,@funcHandleRHSNb,p);
     g4NbEd = integrateVectorised(Nb,curLvl,degreeNeumann,@funcHandleRHSNbVectorised,p);
     n4NbElem = dofU4e(e4ed(NbEd),:);
     neumann = accumarray(n4NbElem(:),g4NbEd(:),[nrNodes+2*nrEdges+nrElems,1]);     
     b = b + neumann;
end
 
%% Include Dirichlet conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = zeros(nrNodes+2*nrEdges+nrElems,1);
firstPoint4ed = 1/3*( c4n(Db(:,2),:) - c4n(Db(:,1),:) ) + c4n(Db(:,1),:);
secondPoint4ed = 2/3*( c4n(Db(:,2),:) - c4n(Db(:,1),:) )+ c4n(Db(:,1),:);
u(fixedNodes) = u_D([c4n(unique(Db),:);firstPoint4ed;secondPoint4ed],p);
b = b - A*u;

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).A = A;
p.level(end).b = b;
p.level(end).x = u;