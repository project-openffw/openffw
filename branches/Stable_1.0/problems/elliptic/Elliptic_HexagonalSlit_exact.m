function p = Elliptic_HexagonalSlit_exact(p)

% PDE definition
p.problem.geom = 'HexagonalSlit';
p.problem.f = @f;
p.problem.g = [];
p.problem.u_D = @u_D;
p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;
p.problem.kappa = @kappa;
p.problem.lambda = @lambda;
p.problem.mu = @mu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Volume force
function z = f(x,y,p)
z = zeros(length(x),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dirichlet boundary values
function z = u_D(x,y,p)
[phi,r] = cart2pol(x,y);
phi( find(phi<0) ) = phi( find(phi<0) ) + 2*pi;
z =  r.^(1/4).*sin(1/4*phi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% elliptic PDE coefficent kappa ( div(kappa*grad_u) )
function z = kappa(x,y,p)

nrElems = length(x);
z = zeros(2,2,nrElems);

for curElem = 1:nrElems 
    z(:,:,curElem) = [1 0;
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
z = r.^(1/4).*sin(1/4*phi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (approximated) gradU
function z = gradU_exact(x,y,p)

h = 1e-10;

z = zeros(length(x),2);
U = u_exact(x,y,p);
Udx = u_exact(x+h,y,p);
Udy = u_exact(x,y+h,p);

z(:,1) = diff([U,Udx],1,2)/h;
z(:,2) = diff([U,Udy],1,2)/h;
