function p = TwoWell_TwoSquare(p)
%author: David Guenther
%Yosida regularization  
%given W*_epsilon by regularisation of W*
%given W_epsilon by conjugation of W*_epsilon

%% PDE definition
p.problem.geom = 'TwoSquare';

p.problem.epsilon = 1e-3;

p = TWgetRegularConj(p);
p = TWgenericNonLinear(p);
%p = TWgetNonLinearRegularConj(p);

p.problem.u_D = @u_D;
%p.problem.f = @f;
%p.problem.f = @f2;
p.problem.f = @f3;




p.problem.hessianU_exact = @hessianU_exact;
p.problem.sigma0 = @sigma0;

syms x y parammu1 parammu2 paramlambda paramt1 real

% Specification of exact solution and differntial Operator
u = x.^6*y.^6*(x-1).^6*(y-1).^6;
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








%% Dirichlet data
function val = u_D(x,y,p)

val = p.problem.u_exact(x,y,p);

%% Volume force
function val = f(x,y,curElem,lvl,p)

val = ones(length(x),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = f2(x,y,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
evalGrad = gradU_exact(x,y,p);
absGrad = (evalGrad(:,1).^2 + evalGrad(:,2).^2).^(1/2);

val = zeros(length(x),1);
t1 = p.problem.t1;
t2 = p.problem.t2;
mu1 = p.problem.mu1;
mu2 = p.problem.mu2;

for k = 1:length(x)
    term1 = ( sin(x(k)).^2 - cos(x(k)).^2 );
    term2 = -( sin(y(k)).^2 - cos(y(k)).^2 );
    term3 = cos(x(k)).^2 * cos(y(k)).^2 + sin(x(k)).^2 * sin(y(k)).^2;
        
    if absGrad(k) <= t1
        val(k) = 0.5 * mu2 * term1 * term2;
    elseif(t1<absGrad(k) && absGrad(k)<t2)
        val(k) = t1 * mu1 * 0.5 * term3^(-0.5) * term1 * term2 - ...
                 t1 * mu2 * 0.25 * term3^(-1.5) * ...
                 sin(x(k))^2 * cos(x(k))^2 * term2^2;
    else    
        val(k) = 0.5 * mu1 * term1 * term2;
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = f3(x,y,curElem,lvl,p)

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


if absFunc ~= 0
DxDW_1 = evalNonLinear2.*varX.*evalFunc(:,1)./absFunc.^2 + ...
    evalNonLinear1.*( squeeze(evalHessian(1,1,:)) - varX.*evalFunc(:,1)./absFunc.^2 );

DyDW_2 = evalNonLinear2.*varY.*evalFunc(:,2)./absFunc.^2 + ...
    evalNonLinear1.*( squeeze(evalHessian(2,2,:)) - varY.*evalFunc(:,2)./absFunc.^2 );

z = -(DxDW_1 + DyDW_2);
else
z=zeros(size(x));
end

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
absFunc = (evalFunc(:,1).^2 + evalFunc(:,2).^2).^(1/2); %|grad(u)|
evalNonLinear = nonLinear(absFunc,curElem,lvl,p);  %DW(|grad(u)|) / |grad(u)|

z = (evalNonLinear*[1,1]).*evalFunc; %DW(|grad(u)|) * (grad(u)*[1,1])/|grad(u)|