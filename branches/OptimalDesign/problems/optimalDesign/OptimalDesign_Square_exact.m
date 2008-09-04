function p = OptimalDesign_Square_exact(p)
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
% p.problem.f = @f2;
p.problem.hessianU_exact = @hessianU_exact;
p.problem.sigma0 = @sigma0;

syms x y parammu1 parammu2 paramlambda paramt1 real

% Specification of exact solution and differntial Operator
u = x*y*(1-x)*(1-y);
gradU = [diff(u,x); diff(u,y)];
hessianU = simple([diff(gradU(1),x); diff(gradU(2),x); diff(gradU(1),y); diff(gradU(2),y)]);

%% exact solution
charU = Matlab4Maple(u);
exec = ['p.problem.u_exact = @(x,y,p)(',charU,');'];
eval(exec,'disp(''err'')');

%% gradU
charGradU = Matlab4Maple(gradU');
exec = ['p.problem.gradU_exact = @(x,y,p)(',charGradU,');'];
eval(exec,'disp(''err'')');

%%  hessianU
charHessianU = Matlab4Maple(hessianU');
exec = ['p.problem.hessianU_exactDummy = @(x,y,p)(',charHessianU,');'];
eval(exec,'disp(''err'')');

%% Volume force
function z = f(x,y,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
evalGrad = gradU_exact(x,y,p);
absGrad = (evalGrad(:,1).^2 + evalGrad(:,2).^2).^(1/2);

z = zeros(length(x),1);
t1 = p.problem.t1;
t2 = p.problem.t2;

for k = 1:length(x)
    if absGrad(k) <= t1
        z(k) = -2.*p.problem.mu2.*(-y(k)+y(k).^2-x(k)+x(k).^2);
    elseif(t1<absGrad(k) && absGrad(k)<t2)
        z(k) = 2.*p.problem.mu2.*p.problem.t1.*x(k).*y(k).*(-1+x(k)).*(-1+y(k)).*(8.*x(k).^2.*y(k).^2-8.*x(k).^2.*y(k)+3.*x(k).^2-8.*x(k).*y(k).^2+8.*x(k).*y(k)-3.*x(k)+3.*y(k).^2-3.*y(k)+1)./(y(k).^2-2.*y(k).^3-4.*x(k).*y(k).^2+8.*x(k).*y(k).^3+y(k).^4-4.*y(k).^4.*x(k)+8.*x(k).^2.*y(k).^2-8.*x(k).^2.*y(k).^3+4.*x(k).^2.*y(k).^4+x(k).^2-4.*x(k).^2.*y(k)-2.*x(k).^3+8.*x(k).^3.*y(k)-8.*x(k).^3.*y(k).^2+x(k).^4-4.*x(k).^4.*y(k)+4.*x(k).^4.*y(k).^2).^(3./2);
    else    
        z(k) = -2.*p.problem.mu1.*(-y(k)+y(k).^2-x(k)+x(k).^2);
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = f2(x,y,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
nonLinear1 = p.problem.nonLinearExactDer;
nonLinear2 = p.problem.nonLinearExactSecDer;
hessianU_exact = p.problem.hessianU_exact;

evalFunc = gradU_exact(x,y,p);
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
function z = u_D(x,y,p)
z = p.problem.u_exact(x,y,p);

%% post-processing hessian
function z = hessianU_exact(x,y,p)

z = p.problem.hessianU_exactDummy(x,y,p);
z = reshape(z,[],4)';
z = reshape(z,2,2,[]);

%% exact stress
function z = sigma0(x,y,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
nonLinear = p.problem.nonLinearExactDer;

evalFunc = gradU_exact(x,y,p);
absFunc = (evalFunc(:,1).^2 + evalFunc(:,2).^2).^(1/2);
evalNonLinear = nonLinear(absFunc,curElem,lvl,p);

z = (evalNonLinear*[1,1]).*evalFunc;
