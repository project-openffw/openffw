function p = TwoWell_Square(p)
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
%p.problem.u_D = @u_D2;
p.problem.f = @f;
p.problem.f0 = @f2;

%% Dirichlet data
function val = u_D(x,y,p)

val = zeros(length(x),1);

function val = u_D2(x,y,p)

%in1 = find(x == 0);
%val(in1) = -3/281216*13^(1/2)*(2*y(in1)-3).^5-1/507*13^(1/2)*(2*y(in1)-3).^3;
%in2 = find(x == 1);
%val(in2) = 1/507*y(in2).^3*13^(1/2)+2/13*y(in2)*13^(1/2);
%in3 = find(y == 0);
%val(in3) = -729/281216*13^(1/2)*(x(in3)-1).^5-9/169*13^(1/2)*(x(in3)-1).^3;
%in4 = find(y == 1.5);
%val(in4) = 1/4056*13^(1/2)*(3*x(in4)+1).^3+1/13*13^(1/2)*(3*x(in4)+1);

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

val = ones(length(x),1);