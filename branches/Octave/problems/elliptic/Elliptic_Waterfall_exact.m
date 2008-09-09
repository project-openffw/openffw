function p = Elliptic_Waterfall_exact(p)
% parameter dependent example on the unit square
%
%   -laplace(u) = f(x,y) in Omega
%            u  = 0 on the boundary
%
% with known exact solution 
%    u = xy(1-x)(1-y)atan(k(sqrt((x-5/4)^2 + (y+1/4)^2)-1))

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


% PDE definition
p.problem.geom = 'Square';
p.problem.f = @f;
p.problem.g = @g;
p.problem.u_D = @u_D;
p.problem.kappa = @kappa;
p.problem.lambda = @lambda;
p.problem.mu = @mu;

p.PDE.k = 15;

% load exact solution
p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Volume force
function z = f(x,y,p)
z = 2.*y.*(1-y).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-1./4.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./4.*x.*y.*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./16.*x.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(3./2).*(32.*x-40).^2./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)-8.*x.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./32.*x.*y.*(1-x).*(1-y).*p.PDE.k.^3./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).*(32.*x-40).^2./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2).^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1)+2.*x.*(1-x).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-1./4.*x.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./4.*x.*y.*(1-x).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./16.*x.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(3./2).*(32.*y+8).^2./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./32.*x.*y.*(1-x).*(1-y).*p.PDE.k.^3./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).*(32.*y+8).^2./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2).^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1);
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

% exact solution
function z = u_exact(x,y,p)
z = x.*y.*(1-x).*(1-y).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradU
function z = gradU_exact(x,y,p)
z = ([[conj(y.*(1-x).*(1-y).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-y.*(1-y).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1)).*x+1./8.*x.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2))+x-x+y-y,conj(x.*(1-x).*(1-y).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-x.*y.*(1-x).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))+1./8.*x.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2))+x-x+y-y]]);


