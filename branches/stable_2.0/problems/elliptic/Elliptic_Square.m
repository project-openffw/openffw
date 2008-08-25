function p = Elliptic_Square(p)
% example on the unit square
%
%   -laplace(u) = 1 in Omega
%            u  = 0 on the boundary

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


% PDE definition
p.problem.geom = 'Square';
p.problem.f = @f;
p.problem.g = @g;
p.problem.u_D = @u_D;
p.problem.kappa = @kappa;
p.problem.Dkappa = @Dkappa;
p.problem.lambda = @lambda;
p.problem.mu = @mu;
return


%% Volume force %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = f(pts,p)
nrPts = size(pts,1);
z = ones(nrPts,1);

%% Dirichlet boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = u_D(pts,p)
nrPts = size(pts,1);
z = zeros(nrPts,1);

%% Neumann boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = g(pts,normals,p)
nrPts = size(pts,1);
z = zeros(nrPts,1);

%% elliptic PDE coefficent kappa ( div(kappa*grad_u) ) %%%%%%%%%%%%%%%%%%%%
function z = kappa(pts,p)
nrPts = size(pts,1);
dim = size(pts,2);
nrPoints = size(pts,1);
z = zeros(dim,dim,nrPts);
for curPoint = 1:nrPoints 
    z(:,:,curPoint) = eye(dim);
end

function z = Dkappa(pts,p)
nrElems = size(pts,1);
% 3rd dimension is derivative with respect to space-dimensions
z = zeros(2,2,2,nrElems);

%% elliptic PDE coefficent lambda ( lambda*grad_u ) %%%%%%%%%%%%%%%%%%%%%%%
function z = lambda(pts,p)
nrPts = size(pts,1);
dim = size(pts,2);
z = zeros(nrPts,dim);
 
%% elliptic PDE coefficent mu ( mu*u ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = mu(pts,p)
nrPts = size(pts,1);
z = zeros(nrPts,1);