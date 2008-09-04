function p = Elasticity_Exact_Template(p)

% PDE definition
p.problem.geom = 'Square';
p.problem.g = @g;
p.problem.u_D = @u_D;
p.problem.sigma_exact = @sigma_exact;

p.PDE.E = 100000;
p.PDE.nu = 0.3;

syms x y paramE paramnu real

u1 = cos((x+1)*(y+1)^2)/1e5;
u2 = sin((x+1))*cos(y+1)/1e5;
u = [u1, u2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exact sigma_exact = C*eps(u)
Du = [diff(u1,x), diff(u1,y);
    diff(u2,x), diff(u2,y)];
epsU = (Du + Du')/2;

paramMu = paramE/(2*(1+paramnu));
paramLambda = paramE * paramnu /( (1+paramnu)*(1-2*paramnu) );
C = paramMu*[2,0,0;0,2,0;0,0,2] + paramLambda*[1,1,0;1,1,0;0,0,0];
sigma_exact = C*[epsU(1,1);epsU(2,2);epsU(2,1)];
sigma_exact = [sigma_exact(1),sigma_exact(3);
                sigma_exact(3),sigma_exact(2)];

charsigma_exact = Matlab4Maple(sigma_exact);
exec = ['p.problem.sigma_exact_dummy = @(x,y,p)(',charsigma_exact,');'];
eval(exec,'disp(''errsigma_exact'')');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Volume force

f = -[diff(sigma_exact(1,1),x) + diff(sigma_exact(1,2),y),...
    diff(sigma_exact(2,1),x) + diff(sigma_exact(2,2),y)];

charF = Matlab4Maple(f);

exec = ['p.problem.f = @(x,y,p)(',charF,');'];
eval(exec,'disp(''errF'')');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exact solution

charU = Matlab4Maple(u);

exec = ['p.problem.u_exact = @(x,y,p)(',charU,');'];
eval(exec,'disp(''errU'')');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dirichlet boundary values
function z = u_D(x,y,p)

z = p.problem.u_exact(x,y,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact Sigma wrapper. (is needed because of Maple -> MATLAB problems when
% x is a vector
function z = sigma_exact(x,y,p)

z = p.problem.sigma_exact_dummy(x,y,p);
z = reshape(z,[],4)';
z = reshape(z,2,2,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Neumann boundary values
function z = g(x,y,n,p)

sigma = p.problem.sigma_exact(x,y,p);
n = reshape(n',2,1,size(n,1));
z = matMul(sigma,n);
z = squeeze(z)';