function p = Elliptic_Square_exact(p)
% example on unit square
%
%   -laplace(u) = 2((x-x^2)+(y-y^2)) in Omega
%            u  = 0 on boundary of Gamma
%
% with known exact solution:
%     u = x(1-x)y(1-y)

% Copyright 2007 Andreas Byfut
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
p.problem.geom = 'Square';%'SquareArbitrary1';%
p.problem.f = @f;
p.problem.g = [];
p.problem.u_D = @u_D;
p.problem.kappa = @kappa;
p.problem.Dkappa = @Dkappa;
p.problem.lambda = @lambda;
p.problem.mu = @mu;

% load exact solution
p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;
return


%% Volume force %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = f(pts,p)
x = pts(:,1);
y = pts(:,2);
z = 2*((x-x.^2)+(y-y.^2));

%% Dirichlet boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = u_D(pts,p)
z = u_exact(pts,p);

%% Neumann boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not given

%% elliptic PDE coefficent kappa ( div(kappa*grad_u) ) %%%%%%%%%%%%%%%%%%%%
function z = kappa(pts,p)
nrPts = size(pts,1);
dim = size(pts,2);
z = zeros(dim,dim,nrPts);
for curPt = 1:nrPts 
    z(:,:,curPt) = eye(dim);
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
 
%% exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = u_exact(pts,p)
x = pts(:,1);
y = pts(:,2);
z = x.*(1-x).*y.*(1-y);

%% gradient of exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = gradU_exact(pts,p)
nrPts = size(pts,1);
x = pts(:,1);
y = pts(:,2);
z = zeros(nrPts,2);
z(:,1) = (1-2*x).*y.*(1-y);
z(:,2) = (1-2*y).*x.*(1-x);