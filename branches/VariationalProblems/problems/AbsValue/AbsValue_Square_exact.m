function p = AbsValue_Square(p)
%author: Lena Noack

%% PDE definition
%p.problem.geom = 'SquareWithoutZero';
p.problem.geom = 'Square';

p.problem.epsilon = 1e-3;

p = AVgetRegularConj(p);
p = AVgenericNonLinear(p);
%p = AVgetNonLinearRegularConj(p);

p.problem.u_D = @u_D;
%p.problem.f = @f;
%p.problem.f = @f2;
p.problem.f = @f3;



%p.problem.u_exact = @u_exact;
%p.problem.gradU_exact = @gradU_exact;
%p.problem.hessianU_exactDummy=@hessianU_exactDummy;
p.problem.sigma0 = @sigma0;
p.problem.hessianU_exact=@hessianU_exact;

syms x y

% Specification of exact solution and differntial operator
%u = x.^2;
u = x*y*(x-1)*(y-1);
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

%function val = u_exact(x,y,p)
%val = sin(x*pi).*sin(y*pi); %min. E(u) = int W(grad u)dx
%
%function val = gradU_exact(x,y,p)
%val = pi*[(cos(x*pi).*sin(y*pi))' (sin(x*pi).*cos(y*pi))'];
%
%function val = hessianU_exactDummy(x,y,p)
%val = pi^2*[(-sin(x*pi).*sin(y*pi))'; (cos(x*pi).*cos(y*pi))'; ...
%            (cos(x*pi).*cos(y*pi))'; (-sin(x*pi).*sin(y*pi))'];

%% Dirichlet data
function val = u_D(x,y,p)

%val = zeros(length(x),1);
val = p.problem.u_exact(x,y,p);

%% Volume force
function val = f(x,y,curElem,lvl,p)

val = ones(length(x),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = f2(x,y,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
evalGrad = gradU_exact(x,y,p);
absGrad = (evalGrad(:,1).^2 + evalGrad(:,2).^2).^(1/2);
nonLinear1 = p.problem.nonLinearExactDer;
nonLinear2 = p.problem.nonLinearExactSecDer;
evalNonLinear1 = nonLinear1(absGrad,curElem,lvl,p); % DW(|Du|)
evalNonLinear2 = nonLinear2(absGrad,curElem,lvl,p); % D2W(|Du|)

val = zeros(length(x),1);

%%D^2W(|Du|)[Du,Du] = W''(|Du|)*(Du*Du)*(Du*Du) + (W'(|Du|)/|Du|)*(Du*Du)

DuDu = evalGrad(:,1).^2 + evalGrad(:,2).^2;

val = evalNonLinear2.*DuDu.*DuDu + evalNonLinear1.*DuDu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = f3(x,y,curElem,lvl,p)
%% D2W(|X|)[Y,Z] = W'' *(Y*X)(Z*X) + (W' /|X|)*(Y*Z)

u_exact = p.problem.u_exact;
gradU_exact = p.problem.gradU_exact;
nonLinear1 = p.problem.nonLinearExactDer;
nonLinear2 = p.problem.nonLinearExactSecDer;
hessianU_exact = p.problem.hessianU_exact;

evalFunc = gradU_exact(x,y,p);                      % Du
absFunc = (evalFunc(:,1).^2 + evalFunc(:,2).^2).^(1/2);
evalNonLinear1 = nonLinear1(absFunc,curElem,lvl,p); % DW(|Du|)
evalNonLinear2 = nonLinear2(absFunc,curElem,lvl,p); % D2W(|Du|)
evalHessian = hessianU_exact(x,y,p);                % D2(u)

varX = evalFunc(:,1).*squeeze(evalHessian(1,1,:)) +...  % 
    evalFunc(:,2).*squeeze(evalHessian(2,1,:));

varY = evalFunc(:,1).*squeeze(evalHessian(1,2,:)) +...
    evalFunc(:,2).*squeeze(evalHessian(2,2,:));


if absFunc ~= 0
%DxDW_1 = evalNonLinear2.*varX.*evalFunc(:,1)./absFunc.^2 + ...
%    evalNonLinear1.*( squeeze(evalHessian(1,1,:)));

%DyDW_2 = evalNonLinear2.*varY.*evalFunc(:,2)./absFunc.^2 + ...
%    evalNonLinear1.*( squeeze(evalHessian(2,2,:)));

                         %-f = W'' * (Du'/|Du|) * D2u * (Du/|Du|)
                         %   + (W' /|Du|) ( Lapl u )
                         
DxDW_1 = evalNonLinear2.*varX.*evalFunc(:,1) + ...
    evalNonLinear1.*( squeeze(evalHessian(1,1,:)));

DyDW_2 = evalNonLinear2.*varY.*evalFunc(:,2) + ...
    evalNonLinear1.*( squeeze(evalHessian(2,2,:)));

z = -(DxDW_1 + DyDW_2);  %-f = W'' * Du' * D2u * Du + (W' /|Du|) ( Lapl u )
            %% D2W(|X|)[Y,Z] = W'' *(Y*X)(Z*X)      + (W' /|X| )*(Y*Z)
else
z=zeros(size(x));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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