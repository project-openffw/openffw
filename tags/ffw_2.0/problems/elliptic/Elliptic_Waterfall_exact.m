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


%% PDE definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.problem.geom = 'Square';%'SquareArbitrary3'; %
p.problem.f = @f;
p.problem.g = @g;
p.problem.u_D = @u_D;
p.problem.kappa = @kappa;
p.problem.Dkappa = @Dkappa;
p.problem.lambda = @lambda;
p.problem.mu = @mu;

% load exact solution
p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;

% paramter for steepness of the gradient of the exact solution
p.PDE.K = 15;

%% specification of exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms x y paramKappa paramDkappa1 paramDkappa2 paramLambda paramMu paramK real 
% u = x*y*(1-x)*(1-y)*atan(paramK*(sqrt((x-5/4)^2 + (y+1/4)^2)-1));

%% elliptic PDE coefficent kappa ( div(kappa*grad_u) ) %%%%%%%%%%%%%%%%%%%%
% paramKappa = [1 0; 0 1];
% paramDkappa1 = [0 0; 0 0];
% paramDkappa2 = [0 0; 0 0];

% Coefficients
% charKappa = Matlab4Maple(paramKappa);
charKappa = '([[1+x-x+y-y,0+x-x+y-y];[0+x-x+y-y,1+x-x+y-y]])';
exec = ['p.problem.kappa_dummy = @(x,y,p)(',charKappa,');'];
eval(exec,'disp(''error initializing Kappa'')');

% charDKappa1 = Matlab4Maple(paramDkappa1);
charDKappa1 = '([[0+x-x+y-y,0+x-x+y-y];[0+x-x+y-y,0+x-x+y-y]])';
exec = ['p.problem.Dkappa1_dummy = @(x,y,p)(',charDKappa1,');'];
eval(exec,'disp(''error initializing DKappa1'')');

% charDKappa2 = Matlab4Maple(paramDkappa2);
charDKappa2 = '([[0+x-x+y-y,0+x-x+y-y];[0+x-x+y-y,0+x-x+y-y]])';
exec = ['p.problem.Dkappa2_dummy = @(x,y,p)(',charDKappa2,');'];
eval(exec,'disp(''error initializing DKappa2'')');

%% elliptic PDE coefficent lambda ( lambda*grad_u ) %%%%%%%%%%%%%%%%%%%%%%%
% paramLambda = [0, 0];
% charLambda = Matlab4Maple(paramLambda);
charLambda = '([[0+x-x+y-y,0+x-x+y-y]])';
exec = ['p.problem.lambda_dummy = @(x,y,p)(',charLambda,');'];
eval(exec,'disp(''error initializing Lambda'')');

%% elliptic PDE coefficent mu ( mu*u ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paramMu = 0;
% charMu = Matlab4Maple(paramMu);
% charMu = [charMu, '+x-x+y-y'];
charMu = '0+x-x+y-y';
exec = ['p.problem.mu_dummy = @(x,y,p)(',charMu,');'];
eval(exec,'disp(''error initializing Mu'')');

%% volume force %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hardcoded right hand side - if symbolic toolbox is not available
charF = '2.*y.*(1-y).*atan(p.PDE.K.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-1./4.*y.*(1-x).*(1-y).*p.PDE.K./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+p.PDE.K.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./4.*x.*y.*(1-y).*p.PDE.K./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+p.PDE.K.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./16.*x.*y.*(1-x).*(1-y).*p.PDE.K./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(3./2).*(32.*x-40).^2./(1+p.PDE.K.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)-8.*x.*y.*(1-x).*(1-y).*p.PDE.K./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)./(1+p.PDE.K.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./32.*x.*y.*(1-x).*(1-y).*p.PDE.K.^3./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).*(32.*x-40).^2./(1+p.PDE.K.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2).^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1)+2.*x.*(1-x).*atan(p.PDE.K.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-1./4.*x.*(1-x).*(1-y).*p.PDE.K./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+p.PDE.K.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./4.*x.*y.*(1-x).*p.PDE.K./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+p.PDE.K.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./16.*x.*y.*(1-x).*(1-y).*p.PDE.K./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(3./2).*(32.*y+8).^2./(1+p.PDE.K.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./32.*x.*y.*(1-x).*(1-y).*p.PDE.K.^3./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).*(32.*y+8).^2./(1+p.PDE.K.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2).^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1)';
% symbolic right hand side - if symbolic toolbox is available
% gradU = [diff(u,x); diff(u,y)];
% KappaGradU = paramKappa * gradU;
% minusDivKappaGradU = -simple(diff(KappaGradU(1)) + diff(KappaGradU(2),y));
% lambdaGradU = paramLambda * gradU;
% muU = paramMu * u;
% paramF = minusDivKappaGradU + lambdaGradU + muU;
% charF = Matlab4Maple(paramF);
exec = ['p.problem.f_dummy = @(x,y,p)(',charF,');'];
eval(exec,'disp(''err'')');

%% exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% charU = Matlab4Maple(u);
charU = 'x.*y.*(1-x).*(1-y).*atan(p.PDE.K.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))';
exec = ['p.problem.u_exact_dummy = @(x,y,p)(',charU,');'];
eval(exec,'disp(''err'')');

%% gradient of exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hardcoded gradient - if symbolic toolbox is not available
charGradU = '([[conj(y.*(1-x).*(1-y).*atan(p.PDE.K.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-y.*(1-y).*atan(p.PDE.K.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1)).*x+1./8.*x.*y.*(1-x).*(1-y).*p.PDE.K./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+p.PDE.K.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2))+x-x+y-y,conj(x.*(1-x).*(1-y).*atan(p.PDE.K.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-x.*y.*(1-x).*atan(p.PDE.K.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))+1./8.*x.*y.*(1-x).*(1-y).*p.PDE.K./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+p.PDE.K.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2))+x-x+y-y]])';
% symbolic right hand side - if symbolic toolbox is available
% charGradU = Matlab4Maple(gradU');
exec = ['p.problem.gradU_exact_dummy = @(x,y,p)(',charGradU,');'];
eval(exec,'disp(''err'')');
return



%% Dirichlet boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = u_D(pts,p)
z = p.problem.u_exact(pts,p);

%% Neumann boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = g(pts,normals,p)
sigma = p.problem.sigma_exact(pts(:,1),pts(:,2),p);
normals = reshape(normals',2,1,size(normals,1));
z = matMul(sigma,normals);
z = squeeze(z)';

%% Wrapper for function handles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = kappa(pts,p)
z = p.problem.kappa_dummy(pts(:,1),pts(:,2),p);
z = reshape(z,[],4)';
z = reshape(z,2,2,[]);

function z = Dkappa(pts,p)
t = p.problem.Dkappa1_dummy(pts(:,1),pts(:,2),p);
t = reshape(t,[],4)';
z(:,:,1,:) = reshape(t,2,2,[]);
t = p.problem.Dkappa2_dummy(pts(:,1),pts(:,2),p);
t = reshape(t,[],4)';
z(:,:,2,:) = reshape(t,2,2,[]);

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
