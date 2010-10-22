function p = TwoWell_Lshape_exact(p)
%author: Lena Noack
%given W** by relaxation of W

%% PDE definition
p.problem.geom = 'Lshape';

p.problem.epsilon = 1e-3;

p = TWgenericNonLinear(p);
p = TWgetNonLinearReg(p);

p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;
p.problem.hessianU_exact = @hessianU_exact;
p.problem.u_D = @u_D;
%p.problem.f = @f;
p.problem.f = @f2;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% exact stress
%% sigma0 = DW(Du)[1,1] 
%% DW(F)[G] = W'1(F)(F-F1,G) + W'2(F)(F-F2,G)
%% DW**(F)[G] = W'**1(F)(F,G) + W'**2(F)(F2,F)(F2,G)
function z = sigma0(x,y,curElem,lvl,p)

F1 = p.problem.F1;
F2 = p.problem.F2;

gradU_exact = p.problem.gradU_exact;

evalFunc = gradU_exact(x,y,p);
gradU1 = evalFunc(:,1); %F(1)
gradU2 = evalFunc(:,2); %F(2)

CONV = p.params.CONV;
if strcmp(CONV,'c')
nonLinearDer1 = p.problem.nonLinearRegDerA;
nonLinearDer2 = p.problem.nonLinearRegDerB;

evalNonLinearDer1 = nonLinearDer1(gradU1,gradU2,curElem,lvl,p);  %W**'1(F)
evalNonLinearDer2 = nonLinearDer2(gradU1,gradU2,curElem,lvl,p);  %W**'2(F)

term1 = (evalNonLinearDer1*[1,1]).*[gradU1 gradU2]; %W'**1(F)(F*[1,1])
term2 = (evalNonLinearDer2*[1,1]).*((F2(1)*gradU1 + F2(2)*gradU2)*[F2(1) F2(2)]); %W'**2(F)(F2*F)(F2*[1,1])

else
nonLinearDer1 = p.problem.nonLinearExactDerA;
nonLinearDer2 = p.problem.nonLinearExactDerB;

evalNonLinearDer1 = nonLinearDer1(gradU1,gradU2,curElem,lvl,p);  %W'1(F)
evalNonLinearDer2 = nonLinearDer2(gradU1,gradU2,curElem,lvl,p);  %W'2(F)

term1 = (evalNonLinearDer1*[1,1]).*[gradU1-F1(1) gradU2-F1(2)]; %W'1(F)*((F-F1)*[1,1])
term2 = (evalNonLinearDer2*[1,1]).*[gradU1-F2(1) gradU2-F2(2)]; %W'2(F)*((F-F2)*[1,1])
end

z= term1 + term2;


%% Dirichlet boundary values
function z = u_D(x,y,p)
z = p.problem.u_exact(x,y,p);

%% Volume force f
function val = f(x,y,curElem,lvl,p)

val = 100*ones(length(x),1);

function z = f2(x,y,curElem,lvl,p)

F1 = p.problem.F1;
F2 = p.problem.F2;
gradU_exact = p.problem.gradU_exact;
evalFunc = gradU_exact(x,y,p);                          %Du
absFunc = (evalFunc(:,1).^2 + evalFunc(:,2).^2).^(1/2); %|Du|
gradU1 = evalFunc(:,1); %F(1)
gradU2 = evalFunc(:,2); %F(2)

CONV = p.params.CONV;
if strcmp(CONV,'c')
nonLinearSecDer1 = p.problem.nonLinearRegSecDerA;
nonLinearSecDer2 = p.problem.nonLinearRegSecDerB;
nonLinearSecDer3 = p.problem.nonLinearRegSecDerC;
evalNonLinear1 = nonLinearSecDer1(gradU1,gradU2,curElem,lvl,p);     %W''**1(|Du|)
evalNonLinear2 = nonLinearSecDer2(gradU1,gradU2,curElem,lvl,p);     %W''**2(|Du|)
evalNonLinear3 = nonLinearSecDer3(gradU1,gradU2,curElem,lvl,p);     %W''**3(|Du|)
else
nonLinearSecDer1 = p.problem.nonLinearExactSecDerA;
nonLinearSecDer2 = p.problem.nonLinearExactSecDerB;
evalNonLinear1 = nonLinearSecDer1(gradU1,gradU2,curElem,lvl,p);     %W''1(|Du|)
evalNonLinear2 = nonLinearSecDer2(gradU1,gradU2,curElem,lvl,p);     %W''2(|Du|)
end

%OD:           D2W(|Du|)*Du*D2u*Du/|Du|^2 + DW(|Du|)/|Du|*(Laplace(u)-Du*D2u*Du/|Du|^2)
%OD:D2W(|X|) = W''(|X|)/|X|^2*(X*Y)(X*Z)  + W'(|X|)/|X|*( Y*Z - (X*Y)(X*Z)/|X|^2 )
%TW(c) : W''**1(|X|)(Laplace(u)) + W''**2(|X|)*Du*D2u*Du + W''**3(|X|)(F2*Z)(F2*Y)
%TW(nc): W''1(|X|)(Laplace(u)) + W''2(|X|)(((X-F1*Y))*((X-F2*Z)) + ((X-F1*Z))*((X-F2*Y)))
%TW(c) : D2W**(|X|)[Y,Z] = W''**1(|X|)(Y*Z) + W''**2(|X|)(X*Y)(X*Z) + W''**3(|X|)(F2*Z)(F2*Y)
%TW(nc): D2W(|X|)[Y,Z]   = W''1(|X|)(Y*Z) + W''2(|X|)(((X-F1*Y))*((X-F2*Z)) + ((X-F1*Z))*((X-F2*Y)))

%(X*Y)(X*Z) = Du*D2u*Du
%Y*Z = Laplace(u)
% wie weiter? Was ist mit (F2*Z),(F2*Y), etc? -> Versuch

%TW(c) : D2W**(|Du|)[Du,Du] = W''**1(|Du|)(Du*Du) + W''**2(|Du|)(Du*Du)(Du*Du) + W''**3(|Du|)(F2*Du)(F2*Du)
%TW(nc): D2W(|Du|)[Du,Du]   = W''1(|Du|)(Du*Du) + W''2(|Du|)*2*(Du-F1)*Du * (Du-F2)*Du

if strcmp(CONV,'c')
    F2Du = F2(1)*evalFunc(:,1)+F2(2)*evalFunc(:,2);
    DuDu = evalFunc(:,1).^2 + evalFunc(:,2).^2;
    D2W = evalNonLinear1.*DuDu + evalNonLinear2.*DuDu.*DuDu + evalNonLinear3.*F2Du.*F2Du;
else
    DuDu = evalFunc(:,1).^2 + evalFunc(:,2).^2;
    DF1D = (evalFunc(:,1)-F1(1)).*evalFunc(:,1) + (evalFunc(:,2)-F1(2)).*evalFunc(:,2);
    DF2D = (evalFunc(:,1)-F2(1)).*evalFunc(:,1) + (evalFunc(:,2)-F2(2)).*evalFunc(:,2);
    D2W = evalNonLinear1.*DuDu + 2*evalNonLinear2.*( DF1D.*DF2D );
end

z = -D2W; 
