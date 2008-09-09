function p = Elliptic_Template(p)
% Template for problem definition 
% with given right hand side.

% Copyright 2007 Jan Reininghaus
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
p.problem.geom = 'Lshape';
p.problem.f = @f;
p.problem.g = @g;
p.problem.u_D = @u_D;
p.problem.kappa = @kappa;
p.problem.lambda = @lambda;
p.problem.mu = @mu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Volume force
function z = f(x,y,p)
z = ones(length(x),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dirichlet boundary values
function z = u_D(x,y,p)
z = zeros(length(x),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary values
function z = g(x,y,n,p)
z = zeros(length(x),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% elliptic PDE coefficent kappa ( div(kappa*grad_u) )
function z = kappa(x,y,p)
nrPoints = length(x);
z = zeros(2,2,nrPoints);
for curPoint = 1:nrPoints 
    z(:,:,curPoint) = [1 0;
                      0 1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% elliptic PDE coefficent lambda ( lambda*grad_u )
function z = lambda(x,y,p)
nrPoints = length(x);
z = zeros(nrPoints,2);
for curPoint = 1:nrPoints 
    z(curPoint,:) = [0 , 0];
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% elliptic PDE coefficent mu ( mu*u )
function z = mu(x,y,p)
z = zeros(length(x),1);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


