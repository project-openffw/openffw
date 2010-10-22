function p = AbsValue_Lshape_exact(p)
%author: Lena Noack

%% PDE definition
p.problem.geom = 'Lshape';

p.problem.epsilon = 1e-3;

p = AVgetRegularConj(p);
p = AVgenericNonLinear(p);

p.problem.u_D = @u_D;
p.problem.f = @f;

%% Dirichlet boundary values
function z = u_D(x,y,p)
z = zeros(length(x),1);

%% Volume force
function val = f(x,y,curElem,lvl,p)

val = ones(length(x),1);