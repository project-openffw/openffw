function p = OptimalDesign_Lshape(p)
%author: David Guenther
%Yosida regularization  
%given W*_epsilon by regularisation of W*
%given W_epsilon by conjugation of W*_epsilon

%% PDE definition
p.problem.geom = 'Lshape';

p.problem.epsilon = 1e-3;

p = ODgenericNonLinear(p);
p = ODgetRegularConj(p);
% p = ODgetNonLinearRegularConj(p);

p.problem.u_D = @u_D;
p.problem.f = @f;

%% Dirichlet data
function val = u_D(points,p)

val = zeros(length(points(:,1)),1);

%% Volume force
function val = f(points,curElem,lvl,p)

val = ones(length(points(:,1)),1);
