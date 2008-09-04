function p = OptimalDesign_SquareSlit(p)
%author: David Guenther
%Yosida regularization  
%given W*_epsilon by regularisation of W*
%given W_epsilon by conjugation of W*_epsilon

%% PDE definition
p.problem.geom = 'SquareSlit';

p.problem.epsilon = 1e-3;

p = ODgenericNonLinear(p);
p = ODgetRegularConj(p);
p = ODgetNonLinearRegularConj(p);

p.problem.u_D = @u_D;
p.problem.f = @f;

%% Dirichlet data
function val = u_D(x,y,p)

val = zeros(length(x),1);

%% Volume force
function val = f(x,y,curElem,lvl,p)

val = ones(length(x),1);
