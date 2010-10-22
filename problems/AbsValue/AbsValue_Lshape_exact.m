function p = AbsValue_Lshape_exact(p)
%author: Lena Noack

%% PDE definition
p.problem.geom = 'Lshape';

p.problem.epsilon = 1e-3;

p = AVgetRegularConj(p);
p = AVgenericNonLinear(p);

p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;
p.problem.hessianU_exact = @hessianU_exact;
p.problem.u_D = @u_D;
%p.problem.f = @f;
%p.problem.f = @f2;
p.problem.f = @f3;
%p.problem.f = @f4;
p.problem.sigma0 = @sigma0;

%% u_exact
function z = u_exact(x,y,p)

[phi,r] = cart2pol(x,y);
ind = find(phi<0);  
phi(ind) = phi(ind)+2*pi*ones(size(ind));
z = r.^(2/3) .* sin(2*phi/3);

%% grad u_exact
function val = gradU_exact(x,y,p)
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
function val = hessianU_exact(x,y,p)
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
function z = sigma0(x,y,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
nonLinear = p.problem.nonLinearExactDer;

evalFunc = gradU_exact(x,y,p);
absFunc = (evalFunc(:,1).^2 + evalFunc(:,2).^2).^(1/2);
evalNonLinear = nonLinear(absFunc,curElem,lvl,p);

z = (evalNonLinear*[1,1]).*evalFunc;

%% Dirichlet boundary values
function z = u_D(x,y,p)
z = p.problem.u_exact(x,y,p);

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
%% Volume force f
function z = f4(x,y,curElem,lvl,p)

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
