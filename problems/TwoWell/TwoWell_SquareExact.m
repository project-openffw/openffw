function p = TwoWell_SquareExact(p)
%author: Lena Noack
%given W** by relaxation of W
%given W*** by conjugation of W**

%% PDE definition
p.problem.geom = 'TwoWellSquare'; %[0,1]x[0,1.5]

p.problem.epsilon = 1e-3;

p = TWgenericNonLinear(p);
p = TWgetNonLinearReg(p);
p = TWgetConj(p);
%p = TWgetConj2(p);

p.problem.u_D = @u_D;
%p.problem.u_D = @u_D2; %noch ändern
p.problem.f = @f;
p.problem.f0 = @f2;



p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;
p.problem.hessianU_exact = @hessianU_exact;
p.problem.sigma0 = @sigma0;

function val = u_exact(x,y,p)
t = (1/sqrt(13)) * (3*(x-1) + 2*y);

in1 = find(t <= 0);
val(in1) = -(3/128)*t(in1).^5 - (1/3)*t(in1).^3;
in2 = find(t > 0);
val(in2) = (1/24)*t(in2).^3 + t(in2);

val=val';

function val = gradU_exact(x,y,p)
t = (1/sqrt(13)) * (3*(x-1) + 2*y);

in1 = find(t <= 0);
value(in1) = -(15/128)*t(in1).^4 - t(in1).^2;
in2 = find(t > 0);
value(in2) = (1/8)*t(in2).^2 + ones(length(in2),1);

val = (1/sqrt(13)) * [3*value' 2*value'];

function val = hessianU_exactDummy(x,y,p)
t = (1/sqrt(13)) * (3*(x-1) + 2*y);

in1 = find(t <= 0);
value(in1) = -(15/32)*t(in1).^3 - 2*t(in1);
in2 = find(t > 0);
value(in2) = (1/4)*t(in2);

val = value*[9/13; 6/13; 6/13; 4/13];


%% Dirichlet data
function val = u_D(x,y,p)

val = p.problem.u_exact(x,y,p);

function val = u_D2(x,y,p)

in1 = find(x == 0);
val(in1) = sqrt(13)*(45657 - 92286*y(in1) + 63144*y(in1).^2 - ...
             15472*y(in1).^3 + 720*y(in1).^4 - 96*y(in1).^5)/843648;
in2 = find(x == 1);
val(in2) = 2*y(in2).*((4/13)*y(in2).^2 + 1)/sqrt(13);
in3 = find(y == 0);
val_dummy = -9*sqrt(13)*(x(in3)-1).^3;
val(in3) = val_dummy.*(81*x(in3).^2-162*x(in3)+1745)/281216;
in4 = find(y == 1.5);
val(in4) = 3*x(in4).*((9/13)*x(in4).^2 + 1)/sqrt(13);

%% Volume force
function val = f(x,y,curElem,lvl,p)

val = zeros(length(x),1);

function val = f2(x,y,curElem,lvl,p)
t = (1/sqrt(13)) * (3*(x-1) + 2*y);

val = -(3/128)*t.^5 - (1/3)*t.^3;

%% post-processing hessian
function z = hessianU_exact(x,y,p)

z = p.problem.hessianU_exactDummy(x,y,p);
z = reshape(z,[],4)';
z = reshape(z,2,2,[]);

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


function z = f3(x,y,curElem,lvl,p)

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
