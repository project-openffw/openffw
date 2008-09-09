function p = Elliptic_Lshape_exact(p)

% PDE definition
p.problem.geom = 'LshapeNeumann';
p.problem.f = @f;
p.problem.g = @g;
p.problem.u_D = @u_D;
p.problem.kappa = @kappa;
p.problem.lambda = @lambda;
p.problem.mu = @mu;

% load exact solution
p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Volume force
function z = f(x,y,p)
z = zeros(length(x),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Neumann boundary values (Stress)
function z = g(x,y,normals,p)
[phi,r] = cart2pol(x,y);
phi = 2*pi*(phi<-eps)+phi; % This yields 0 <= phi < 2 pi
for j = 1:length(x)
    rotNormal = [cos(phi(j)) sin(phi(j));-sin(phi(j)) cos(phi(j))]*normals(j,:)';
    z(j) = 2/3*r(j)^(-1/3)*[sin(2/3*phi(j)), cos(2/3*phi(j))]*rotNormal;
end
z = z';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dirichlet boundary values
function z = u_D(x,y,p)
z = u_exact(x,y,p);

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
[phi,r] = cart2pol(x,y);
phi( find(phi<0) ) = phi( find(phi<0) ) + 2*pi;
z =  r.^(2/3).*sin(2/3*phi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (approximated) gradU
function z = gradU_exact(x,y,p)

z = zeros(length(x),2);

[phi,r] = cart2pol(x,y);
phi = 2*pi*(phi<-eps)+phi; % This yields 0 <= phi < 2 pi
for j = 1:length(x)
    rotNormal = [cos(phi(j)) sin(phi(j));-sin(phi(j)) cos(phi(j))];
    z(j,:) = 2/3*r(j)^(-1/3)*[sin(2/3*phi(j)), cos(2/3*phi(j))]*rotNormal;
end


