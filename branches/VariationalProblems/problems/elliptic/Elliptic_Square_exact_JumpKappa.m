function p = Elliptic_Square_exact_JumpKappa(p)

% PDE definition
p.problem.geom = 'SquareBig';
p.problem.g = [];
p.problem.u_D = @u_D;
p.problem.kappa = @kappa;
p.problem.f = @f;

% load exact solution
p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;

syms x y kappa lambda mu real



% Specification of exact solution and differntial Operator
u = sin(x^3)*cos(y^3);;
lambda = [0, 0];
mu = 0;

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

lambdaGradU = lambda*gradU;
muU = mu*u;

% f = minusDivKappaGradU + lambdaGradU + muU;
% f = simple(real(f));
% charF = Matlab4Maple(f);
% exec = ['p.problem.f = @(x,y,p)(',charF,');'];
% eval(exec,'disp(''err'')');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = f(x,y,p)
z = zeros(length(x),1);

I = find(abs(x)<=1);
J = find(abs(x)>1);
% I = 1:length(x);

z(I) = 9.*sin(x(I).^3).*x(I).^4.*cos(y(I).^3)-6.*cos(x(I).^3).*x(I).*cos(y(I).^3)+9.*sin(x(I).^3).*cos(y(I).^3).*y(I).^4+6.*sin(x(I).^3).*sin(y(I).^3).*y(I);
z(J) = 36.*sin(x(J).^3).*x(J).^4.*cos(y(J).^3)-24.*cos(x(J).^3).*x(J).*cos(y(J).^3)+36.*sin(x(J).^3).*cos(y(J).^3).*y(J).^4+24.*sin(x(J).^3).*sin(y(J).^3).*y(J);

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

% z = p.problem.kappa_dummy(x,y,p);
% z = reshape(z,[],4)';
% z = reshape(z,2,2,[]);
z = zeros(2,2,length(x));
for k = 1:length(x)
    if abs(x(k)) <= 1
        z(:,:,k) = [1 0;0 1];
    else
        z(:,:,k) = [2 0;0 2];
    end
end
