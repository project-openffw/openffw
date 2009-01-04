function p = OptimalDesign_Square_exact2(p)
%author: David Guenther
%Yosida regularization  
%given W*_epsilon by regularisation of W*
%given W_epsilon by conjugation of W*_epsilon

%% PDE definition
p.problem.geom = 'Square';

p.problem.epsilon = 1e-3;

p = ODgenericNonLinear(p);
p = ODgetRegularConj(p);
% p = ODgetNonLinearRegularConj(p);

p.problem.u_D = @u_D;
p.problem.f = @f;
p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;
p.problem.hessianU_exact = @hessianU_exact;
p.problem.sigma0 = @sigma0;

syms x y parammu1 parammu2 paramlambda paramt1 real

% Specification of exact solution and differntial Operator
% u = sin(2*pi*x)*sin(2*pi*y);
u = x*y*(1-x)*(1-y)*atan(15*(sqrt((x-5/4)^2 + (y+1/4)^2)-1));

gradU = [diff(u,x); diff(u,y)];
hessianU = simple([diff(gradU(1),x); diff(gradU(2),x); diff(gradU(1),y); diff(gradU(2),y)]);

%% exact solution
charU = Matlab4Maple(u);
exec = ['p.problem.u_exact_dummy = @(x,y,p)(',charU,');'];
eval(exec,'disp(''err'')');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% gradU
charGradU = Matlab4Maple(gradU');
exec = ['p.problem.gradU_exact_dummy = @(x,y,p)(',charGradU,');'];
eval(exec,'disp(''err'')');

%% hessianU
charHessianU = Matlab4Maple(hessianU');
exec = ['p.problem.hessianU_exact_dummy = @(x,y,p)(',charHessianU,');'];
eval(exec,'disp(''err'')');

%% Volume force
function z = f(points,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
nonLinear1 = p.problem.nonLinearExactDer;
nonLinear2 = p.problem.nonLinearExactSecDer;
hessianU_exact = p.problem.hessianU_exact;

evalFunc = gradU_exact(points,p);
absFunc = (evalFunc(:,1).^2 + evalFunc(:,2).^2).^(1/2);
evalNonLinear1 = nonLinear1(absFunc,curElem,lvl,p);
evalNonLinear2 = nonLinear2(absFunc,curElem,lvl,p);
evalHessian = hessianU_exact(x,y,p);

varX = evalFunc(:,1).*squeeze(evalHessian(1,1,:)) +...
    evalFunc(:,2).*squeeze(evalHessian(2,1,:));

varY = evalFunc(:,1).*squeeze(evalHessian(1,2,:)) +...
    evalFunc(:,2).*squeeze(evalHessian(2,2,:));

DxDW_1 = evalNonLinear2.*varX.*evalFunc(:,1)./absFunc.^2 + ...
    evalNonLinear1.*( squeeze(evalHessian(1,1,:)) - varX.*evalFunc(:,1)./absFunc.^2 );

DyDW_2 = evalNonLinear2.*varY.*evalFunc(:,2)./absFunc.^2 + ...
    evalNonLinear1.*( squeeze(evalHessian(2,2,:)) - varY.*evalFunc(:,2)./absFunc.^2 );

z = -(DxDW_1 + DyDW_2);

%% Dirichlet boundary values
function z = u_D(points,p)
z = p.problem.u_exact(points,p);

%% post-processing hessian
function z = hessianU_exact(points,p)

z = p.problem.hessianU_exactDummy(points,p);
z = reshape(z,[],4)';
z = reshape(z,2,2,[]);

%% exact stress
function z = sigma0(points,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
nonLinear = p.problem.nonLinearExactDer;

evalFunc = gradU_exact(points,p);
absFunc = (evalFunc(:,1).^2 + evalFunc(:,2).^2).^(1/2);
evalNonLinear = nonLinear(absFunc,curElem,lvl,p);

z = (evalNonLinear*[1,1]).*evalFunc;

%% post-processing gradient
function z = gradU_exact(points,p)

x = points(:,1);
y = points(:,2);

z = p.problem.gradU_exact_dummy(x,y,p);

%% post-processing displacement
function z = u_exact(points,p)

x = points(:,1);
y = points(:,2);

z = p.problem.u_exact_dummy(x,y,p);
