function p = Elasticity_Template_Exact(p)

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


%% Problem definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDE definition
p.problem.geom = 'Square';
p.problem.f = @f;
p.problem.g = @g;
p.problem.u_D = @u_D;

% load exact solution
p.problem.u_exact = @u_exact;
p.problem.sigma_exact = @sigma_exact;

% set Lamé parameters
p.PDE.E = 100000;
p.PDE.nu = 0.3;

%% specification of exact solution, etc. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x y paramE paramnu real
u1 = cos((x+1)*(y+1)^2)/1e5;
u2 = sin((x+1))*cos(y+1)/1e5;
u = [u1, u2];

% exact sigma_exact = C*eps(u)
Du = [diff(u1,x), diff(u1,y);
      diff(u2,x), diff(u2,y)];
epsU = (Du + Du')/2;
paramMu = paramE/(2*(1+paramnu));
paramLambda = paramE * paramnu /( (1+paramnu)*(1-2*paramnu) );
C = paramMu*[2,0,0;0,2,0;0,0,2] + paramLambda*[1,1,0;1,1,0;0,0,0];
sigma_exact_dummy = C*[epsU(1,1);epsU(2,2);epsU(2,1)];
sigma_exact_dummy = [sigma_exact_dummy(1),sigma_exact_dummy(3);
                     sigma_exact_dummy(3),sigma_exact_dummy(2)];
charsigma_exact = Matlab4Maple(sigma_exact_dummy);
exec = ['p.problem.sigma_exact_dummy = @(x,y,p)(',charsigma_exact,');'];
eval(exec,'disp(''errsigma_exact'')');

% Volume force
f_dummy = -[diff(sigma_exact_dummy(1,1),x) + diff(sigma_exact_dummy(1,2),y),...
            diff(sigma_exact_dummy(2,1),x) + diff(sigma_exact_dummy(2,2),y)];
charF = Matlab4Maple(f_dummy);
exec = ['p.problem.f = @(x,y,p)(',charF,');'];
eval(exec,'disp(''errF'')');

% exact solution
charU = Matlab4Maple(u);
exec = ['p.problem.u_exact = @(x,y,p)(',charU,');'];
eval(exec,'disp(''errU'')');
return



%% Dirichlet boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = u_D(pts,p)
z = p.problem.u_exact(pts,p);

%% Neumann boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = g(pts,n,p)
sigma = p.problem.sigma_exact(pts,p);
n = reshape(n',2,1,size(n,1));
z = matMul(sigma,n);
z = squeeze(z)';



%% Wrapper for function handles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = f(pts,p)
z = p.problem.f_dummy(pts(:,1),pts(:,2),p);

function z = u_exact(pts,p)
z = p.problem.u_exact_dummy(pts(:,1),pts(:,2),p);

function z = sigma_exact(pts,p)
z = p.problem.sigma_exact_dummy(pts(:,1),pts(:,2),p);
z = reshape(z,[],4)';
z = reshape(z,2,2,[]);