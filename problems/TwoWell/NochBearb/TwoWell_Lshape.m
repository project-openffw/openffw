function p = TwoWell_Lshape(p)
%author: Lena Noack
%given W** by relaxation of W

%% PDE definition
p.problem.geom = 'Lshape';

p.problem.epsilon = 1e-3;

p = TWgenericNonLinear(p);
p = TWgetNonLinearReg(p);

p.problem.u_D = @u_D;
p.problem.f = @f;

%% Dirichlet data
function val = u_D(x,y,p)

val = zeros(length(x),1);

%% Volume force
function val = f(x,y,curElem,lvl,p)

val = ones(length(x),1);
