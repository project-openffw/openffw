function p = Elliptic_Waterfall_exact(p)

% PDE definition
p.problem.geom = 'Square';
p.problem.g = [];
p.problem.u_D = @u_D;
p.problem.kappa = @kappa;

% load exact solution
p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;

p.PDE.k = 15;

syms x y kappa lambda mu paramk real 

% Specification of exact solution and differntial Operator
u = x*y*(1-x)*(1-y)*atan(paramk*(sqrt((x-5/4)^2 + (y+1/4)^2)-1));

lambda = [0, 0];
mu = 0;
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
% hardcoded right hand side - we do not assume that the symbolic toolbox is
% always available. 

% gradU = [diff(u,x); diff(u,y)];
% KappaGradU = kappa * gradU;
% minusDivKappaGradU = -simple(diff(KappaGradU(1)) + diff(KappaGradU(2),y));
% 
% lambdaGradU = lambda*gradU;
% muU = mu*u;
% 
% f = minusDivKappaGradU + lambdaGradU + muU;
% charF = Matlab4Maple(f);

charF = '2.*y.*(1-y).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-1./4.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./4.*x.*y.*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./16.*x.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(3./2).*(32.*x-40).^2./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)-8.*x.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./32.*x.*y.*(1-x).*(1-y).*p.PDE.k.^3./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).*(32.*x-40).^2./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2).^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1)+2.*x.*(1-x).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-1./4.*x.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./4.*x.*y.*(1-x).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./16.*x.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(3./2).*(32.*y+8).^2./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./32.*x.*y.*(1-x).*(1-y).*p.PDE.k.^3./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).*(32.*y+8).^2./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2).^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1)';

exec = ['p.problem.f = @(x,y,p)(',charF,');'];
eval(exec,'disp(''err'')');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exact solution
charU = Matlab4Maple(u);
exec = ['p.problem.u_exact = @(x,y,p)(',charU,');'];
eval(exec,'disp(''err'')');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  gradU
% charGradU = Matlab4Maple(gradU');
charGradU = '([[conj(y.*(1-x).*(1-y).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-y.*(1-y).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1)).*x+1./8.*x.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2))+x-x+y-y,conj(x.*(1-x).*(1-y).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-x.*y.*(1-x).*atan(p.PDE.k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))+1./8.*x.*y.*(1-x).*(1-y).*p.PDE.k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+p.PDE.k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2))+x-x+y-y]])';
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


