function p = Elliptic_Template_exact(p)
% Template for problem definition
% for given exact solution u
% using the 'symbolic toolbox'

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


%% PDE definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.problem.geom = 'Square';
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


%% specification of exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x y paramKappa paramDkappa1 paramDkappa2 paramLambda paramMu real
u = sin(x^3)*cos(y^pi)+x^8-y^9+x^6*y^10;

%% elliptic PDE coefficent kappa ( div(kappa*grad_u) ) %%%%%%%%%%%%%%%%%%%%
paramKappa = [1 0; 0 1];
paramDkappa1 = [0 0; 0 0];
paramDkappa2 = [0 0; 0 0];
charKappa = Matlab4Maple(paramKappa);
exec = ['p.problem.kappa_dummy = @(x,y,p)(',charKappa,');'];
eval(exec,'disp(''error initializing Kappa'')');
charDKappa1 = Matlab4Maple(paramDkappa1);
exec = ['p.problem.Dkappa1_dummy = @(x,y,p)(',charDKappa1,');'];
eval(exec,'disp(''error initializing DKappa1'')');
charDKappa2 = Matlab4Maple(paramDkappa2);
exec = ['p.problem.Dkappa2_dummy = @(x,y,p)(',charDKappa2,');'];
eval(exec,'disp(''error initializing DKappa2'')');

%% elliptic PDE coefficent lambda ( lambda*grad_u ) %%%%%%%%%%%%%%%%%%%%%%%
paramLambda = [0, 0];
charLambda = Matlab4Maple(paramLambda);
exec = ['p.problem.lambda_dummy = @(x,y,p)(',charLambda,');'];
eval(exec,'disp(''error initializing Lambda'')');

%% elliptic PDE coefficent mu ( mu*u ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramMu = 0;
charMu = Matlab4Maple(paramMu);
charMu = [charMu, '+x-x+y-y'];
exec = ['p.problem.mu_dummy = @(x,y,p)(',charMu,');'];
eval(exec,'disp(''error initializing Mu'')');

%% volume force %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradU = [diff(u,x); diff(u,y)];
KappaGradU = paramKappa * gradU;
minusDivKappaGradU = -simple(diff(KappaGradU(1),x) + diff(KappaGradU(2),y));
lambdaGradU = paramLambda * gradU;
muU = paramMu * u;
paramF = minusDivKappaGradU + lambdaGradU + muU;
charF = Matlab4Maple(paramF);
exec = ['p.problem.f_dummy = @(x,y,p)(',charF,');'];
eval(exec,'disp(''err'')');

%% exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
charU = Matlab4Maple(u);
exec = ['p.problem.u_exact_dummy = @(x,y,p)(',charU,');'];
eval(exec,'disp(''err'')');

%% gradient of exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
charGradU = Matlab4Maple(gradU');
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