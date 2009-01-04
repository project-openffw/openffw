function p = OptimalDesign_Lshape_exact(p)
%author: David Guenther
%Yosida regularization  
%given W*_epsilon by regularisation of W*
%given W_epsilon by conjugation of W*_epsilon

%% PDE definition
p.problem.geom = 'Lshape';

p.problem.epsilon = 1e-3;

p = ODgenericNonLinear(p);
p = ODgetRegularConj(p);
% p = ODgetNonLinearRegularConj(p);

p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;
p.problem.hessianU_exact = @hessianU_exact;
p.problem.u_D = @u_D;
p.problem.f = @f;
p.problem.sigma0 = @sigma0;

%% u_exact
function z = u_exact(points,p)

x = points(:,1);
y = points(:,2);

[phi,r] = cart2pol(x,y);
ind = find(phi<0);  
phi(ind) = phi(ind)+2*pi*ones(size(ind));
z = r.^(2/3) .* sin(2*phi/3);

%% grad u_exact
function val = gradU_exact(points,p)

x = points(:,1);
y = points(:,2);

[phi,r] = cart2pol(x,y);
ind = find(phi<0);  
phi(ind) = phi(ind)+2*pi*ones(size(ind));

val = zeros(length(r),2);

for k = 1:length(r)
    matrix = [   cos(phi(k)) ,      -sin(phi(k))./r(k);...
                 sin(phi(k)),   cos(phi(k))./r(k)    ];

    vector  = [(2/3) * r(k).^(-1/3) .* ( sin(2*phi(k)/3)),...
               (2/3) * r(k).^(2/3) .* ( cos(2*phi(k)/3))];
    val(k,:) = matrix*vector';
end

%% Hessian of u 
function val = hessianU_exact(points,p)

x = points(:,1);
y = points(:,2);

[phi,r] = cart2pol(x,y);
ind = find(phi<0);  
phi(ind) = phi(ind)+2*pi*ones(size(ind));

dF1dR = 2/9*sin(1/3*phi)./r.^(4/3);
dF1dP = -2/9*cos(1/3*phi)./r.^(1/3);
dF2dR = -2/9*cos(1/3*phi)./r.^(4/3);
dF2dP = -2/9*sin(1/3*phi)./r.^(1/3);

val = zeros(2,2,length(x));

for k = 1:length(x)
    matrix = [   cos(phi(k)) ,      -sin(phi(k))./r(k);...
                 sin(phi(k)) ,   cos(phi(k))./r(k)    ];
    dF = [dF1dR(k), dF2dR(k);
          dF1dP(k), dF2dP(k)  ];
      
    val(:,:,k) = (matrix*dF)';
end

%% sigma_0: p = DW(\gradU) = W'(|\gradU|)/|gradU|*gradU
function z = sigma0(points,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
nonLinear = p.problem.nonLinearExactDer;

x = points(:,1);
y = points(:,2);

evalFunc = gradU_exact(x,y,p);
absFunc = (evalFunc(:,1).^2 + evalFunc(:,2).^2).^(1/2);
evalNonLinear = nonLinear(absFunc,curElem,lvl,p);

z = (evalNonLinear*[1,1]).*evalFunc;

%% Dirichlet boundary values
function z = u_D(points,p)

z = p.problem.u_exact(points,p);

%% Volume force f
function z = f(points,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
nonLinear1 = p.problem.nonLinearExactDer;
nonLinear2 = p.problem.nonLinearExactSecDer;
hessianU_exact = p.problem.hessianU_exact;

evalFunc = gradU_exact(points,p);
absFunc = (evalFunc(:,1).^2 + evalFunc(:,2).^2).^(1/2);
evalNonLinear1 = nonLinear1(absFunc,curElem,lvl,p);
evalNonLinear2 = nonLinear2(absFunc,curElem,lvl,p);
evalHessian = hessianU_exact(points,p);

varX = evalFunc(:,1).*squeeze(evalHessian(1,1,:)) +...
    evalFunc(:,2).*squeeze(evalHessian(2,1,:));

varY = evalFunc(:,1).*squeeze(evalHessian(1,2,:)) +...
    evalFunc(:,2).*squeeze(evalHessian(2,2,:));

DxDW_1 = evalNonLinear2.*varX.*evalFunc(:,1)./absFunc.^2 + ...
    evalNonLinear1.*( squeeze(evalHessian(1,1,:)) - varX.*evalFunc(:,1)./absFunc.^2 );

DyDW_2 = evalNonLinear2.*varY.*evalFunc(:,2)./absFunc.^2 + ...
    evalNonLinear1.*( squeeze(evalHessian(2,2,:)) - varY.*evalFunc(:,2)./absFunc.^2 );

z = -(DxDW_1 + DyDW_2);