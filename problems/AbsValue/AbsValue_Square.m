function p = AbsValue_Square(p)
%author: Lena Noack

%% PDE definition
%p.problem.geom = 'SquareWithoutZero';
p.problem.geom = 'Square';

p.problem.epsilon = 1e-3;

p = AVgenericNonLinear(p);
p = AVgetRegularConj(p);
%p = AVgetNonLinearRegularConj(p);

p.problem.u_D = @u_D;
p.problem.f = @f;


%% Dirichlet data
function val = u_D(x,y,p)

val = zeros(length(x),1);

%% Volume force
function val = f(x,y,curElem,lvl,p)

val = ones(length(x),1);

