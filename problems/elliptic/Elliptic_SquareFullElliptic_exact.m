function p = Elliptic_SquareFullElliptic_exact(p)
% full elliptic example on the unit square
%
%   -laplace(u)+ lambda grad(u) + mu u = f(x,y) in Omega
%            u  = u_D(x,y) on the boundary
%
% with known exact solution 
%     u = sin(x^3)cos(y^pi)+x^8-y^9+x^6y^10 

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


%% PDE definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDE definition
p.problem.geom = 'Square';%'SquareArbitrary2';%
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

%% specification of exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms x y paramKappa paramDkappa1 paramDkappa2 paramLambda paramMu real
% u = sin(x^3)*cos(y^pi)+x^8-y^9+x^6*y^10;

%% elliptic PDE coefficent kappa ( div(kappa*grad_u) ) %%%%%%%%%%%%%%%%%%%%
% paramKappa = [1 0; 0 1];
% paramDkappa1 = [0 0; 0 0];
% paramDkappa2 = [0 0; 0 0];
% charKappa = Matlab4Maple(paramKappa);
charKappa = '([[1+x-x+y-y,0+x-x+y-y];[0+x-x+y-y,1+x-x+y-y]])';
exec = ['p.problem.kappa_dummy = @(x,y,p)(',charKappa,');'];
eval(exec,'disp(''error initializing Kappa'')');
% charDKappa1 = Matlab4Maple(paramDkappa1);
charDKappa1 = '([[0+x-x+y-y,0+x-x+y-y];[0+x-x+y-y,0+x-x+y-y]])';
exec = ['p.problem.Dkappa1_dummy = @(x,y,p)(',charDKappa1,');'];
eval(exec,'disp(''error initializing DKappa1'')');
% charDKappa2 = Matlab4Maple(paramDkappa2);
charDKappa2= '([[0+x-x+y-y,0+x-x+y-y];[0+x-x+y-y,0+x-x+y-y]])';
exec = ['p.problem.Dkappa2_dummy = @(x,y,p)(',charDKappa2,');'];
eval(exec,'disp(''error initializing DKappa2'')');

%% elliptic PDE coefficent lambda ( lambda*grad_u ) %%%%%%%%%%%%%%%%%%%%%%%
% paramLambda = [5*sin(x+y), 6*cos(x+y)];
% charLambda = Matlab4Maple(paramLambda);
charLambda = '([[5.*sin(x+y)+x-x+y-y,6.*cos(x+y)+x-x+y-y]])';
exec = ['p.problem.lambda_dummy = @(x,y,p)(',charLambda,');'];
eval(exec,'disp(''error initializing Lambda'')');

%% elliptic PDE coefficent mu ( mu*u ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paramMu = 7;
% charMu = Matlab4Maple(paramMu);
% charMu = [charMu, '+x-x+y-y'];
charMu = '7+x-x+y-y';
exec = ['p.problem.mu_dummy = @(x,y,p)(',charMu,');'];
eval(exec,'disp(''error initializing Mu'')');

%% volume force %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gradU = [diff(u,x); diff(u,y)];
% KappaGradU = paramKappa * gradU;
% minusDivKappaGradU = -simple(diff(KappaGradU(1)) + diff(KappaGradU(2),y));
% lambdaGradU = paramLambda * gradU;
% muU = paramMu * u;
% paramF = minusDivKappaGradU + lambdaGradU + muU;
% charF = Matlab4Maple(paramF);
charF = '-(56+90.*y.^8).*x.^6-(-9.*sin(x.^3).*cos(y.^pi)+30.*y.^10).*x.^4-6.*cos(x.^3).*x.*cos(y.^pi)+sin(x.^3).*cos(y.^pi).*(y.^pi).^2.*pi.^2./y.^2+sin(x.^3).*sin(y.^pi).*y.^pi.*pi.^2./y.^2-sin(x.^3).*sin(y.^pi).*y.^pi.*pi./y.^2+72.*y.^7+5.*sin(x+y).*(3.*cos(x.^3).*x.^2.*cos(y.^pi)+8.*x.^7+6.*x.^5.*y.^10)+6.*cos(x+y).*(-sin(x.^3).*sin(y.^pi).*y.^pi.*pi./y-9.*y.^8+10.*x.^6.*y.^9)+7.*sin(x.^3).*cos(y.^pi)+7.*x.^8-7.*y.^9+7.*x.^6.*y.^10';
exec = ['p.problem.f_dummy = @(x,y,p)(',charF,');'];
eval(exec,'disp(''err'')');

%% exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% charU = Matlab4Maple(u);
charU = 'sin(x.^3).*cos(y.^pi)+x.^8-y.^9+x.^6.*y.^10';
exec = ['p.problem.u_exact_dummy = @(x,y,p)(',charU,');'];
eval(exec,'disp(''err'')');

%% gradient of exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% charGradU = Matlab4Maple(gradU');
charGradU = '([[8.*x.^7+6.*x.^5.*y.^10+3.*cos(x.^3).*x.^2.*cos(conj(y.^pi))+x-x+y-y,-9.*y.^8+10.*x.^6.*y.^9-sin(x.^3).*pi./y.*conj(sin(y.^pi).*y.^pi)+x-x+y-y]])';
exec = ['p.problem.gradU_exact_dummy = @(x,y,p)(',charGradU,');'];
eval(exec,'disp(''err'')');
return



%% Dirichlet boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = u_D(pts,p)
z = p.problem.u_exact(pts,p);

%% Neumann boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = g(pts,normals,p)
sigma = p.problem.sigma_exact(pts,p);
normals = reshape(normals',2,1,size(normals,1));
z = matMul(sigma,normals);
z = squeeze(z)';



%% Wrapper for function handles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = kappa(pts,p)
z = p.problem.kappa_dummy(pts(:,1),pts(:,2),p);
z = reshape(z,[],4)';
z = reshape(z,2,2,[]);

function z = lambda(pts,p)
z = p.problem.lambda_dummy(pts(:,1),pts(:,2),p);

function z = mu(pts,p)
z = p.problem.mu_dummy(pts(:,1),pts(:,2),p);

function z = f(pts,p)
z = p.problem.f_dummy(pts(:,1),pts(:,2),p);

function z = u_exact(pts,p)
z = p.problem.u_exact_dummy(pts(:,1),pts(:,2),p);

function z = gradU_exact(pts,p)
z = p.problem.gradU_exact_dummy(pts(:,1),pts(:,2),p);

function z = Dkappa(pts,p)
t = p.problem.Dkappa1_dummy(pts(:,1),pts(:,2),p);
t = reshape(t,[],4)';
z(:,:,1,:) = reshape(t,2,2,[]);
t = p.problem.Dkappa2_dummy(pts(:,1),pts(:,2),p);
t = reshape(t,[],4)';
z(:,:,2,:) = reshape(t,2,2,[]);
