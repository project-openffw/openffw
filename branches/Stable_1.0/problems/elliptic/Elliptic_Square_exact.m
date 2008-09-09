function p = Elliptic_Square_exact(p)

% PDE definition
p.problem.geom = 'Square';
p.problem.g = [];
p.problem.u_D = @u_D;
p.problem.kappa = @kappa;

% load exact solution
p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;

syms x y kappa lambda mu real



% Specification of exact solution and differntial Operator
u = sin(x^3)*cos(y^pi)+x^8-y^9+x^6*y^10;
lambda = [5*sin(x+y), 6*cos(x+y)];
mu = 7;
kappa = [1 0; 0 1];



% Coefficients
charKappa = Matlab4Maple(kappa);
exec = ['p.problem.kappa_dummy = @(x,y,p)(',charKappa,');'];
eval(exec,'disp(''error initializing Kappa'')');

charLambda = Matlab4Maple(lambda);
exec = ['p.problem.lambda = @(x,y,p)(',charLambda,');'];
eval(exec,'disp(''error initializing Lambda'')');

charMu = Matlab4Maple(mu);
charMu = [charMu, '+x-x+y-y'];
exec = ['p.problem.mu = @(x,y,p)(',charMu,');'];
eval(exec,'disp(''error initializing Mu'')');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Volume force
gradU = [diff(u,x); diff(u,y)];
KappaGradU = kappa * gradU;
minusDivKappaGradU = -simple(diff(KappaGradU(1)) + diff(KappaGradU(2),y));

lambdaGradU = lambda*gradU;
muU = mu*u;

f = minusDivKappaGradU + lambdaGradU + muU;
f = simple(real(f));
charF = Matlab4Maple(f);
exec = ['p.problem.f = @(x,y,p)(',charF,');'];
eval(exec,'disp(''err'')');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exact solution
charU = Matlab4Maple(u);
exec = ['p.problem.u_exact = @(x,y,p)(',charU,');'];
eval(exec,'disp(''err'')');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  gradU
charGradU = Matlab4Maple(gradU');
exec = ['p.problem.gradU_exact = @(x,y,p)(',charGradU,');'];
eval(exec,'disp(''err'')');

% Dirichlet boundary values
function z = u_D(x,y,p)
z = p.problem.u_exact(x,y,p);

% Neumann boundary values
function z = g(x,y,n,p)
sigma = p.problem.sigma_exact(x,y,p);
n = reshape(n',2,1,size(n,1));
z = matMul(sigma,n);
z = squeeze(z)';

% Exact Sigma wrapper. (is needed because of Maple -> MATLAB problems when
% x is a vector
function z = kappa(x,y,p)

z = p.problem.kappa_dummy(x,y,p);
z = reshape(z,[],4)';
z = reshape(z,2,2,[]);
