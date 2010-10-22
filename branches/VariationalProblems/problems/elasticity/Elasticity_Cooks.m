function p = Elasticity_Cooks(p)

% PDE definition
p.problem.geom = 'Cooks';
p.problem.f = @f;
p.problem.g = @g;
p.problem.u_D = @u_D;
p.problem.u_exact = [];
p.PDE.E = 2900; 
p.PDE.nu = 0.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Volume force
function z = f(x,y,p)
z = zeros(size(x,1),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z = g(x,y,n,p)
z = zeros(size(x,1),2);
z(find(n(:,1)==1),2) = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dirichlet boundary values
function z = u_D(x,y,p)
z = zeros(length(x),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
