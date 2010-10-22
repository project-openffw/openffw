function p = genericNonLinear(p)

p.problem.nonLinearConjExact = @nonLinearConjExact;
p.problem.nonLinearExact = @nonLinearExact;
p.problem.nonLinearExactDer = @nonLinearExactDer;
p.problem.nonLinearExactSecDer = @nonLinearExactSecDer;

%% exact conjugate nonlinear operator W*(x) (not regularised)
function z = nonLinearConjExact(x,curElem,lvl,p)

mu1 = p.problem.mu1;
mu2 = p.problem.mu2;
t1 = p.problem.t1;
t2 = p.problem.t2;

z = zeros(length(x),1);

index1 = find( x <= mu2*t1 );
z(index1) = 1/(2*mu2) * x(index1).^2;

index2 = find( x > mu2*t2 );
z(index2) = 1/(2*mu1) * x(index2).^2 - 1/2 * mu1 * t2^2 + 1/2 * mu2 * t1^2;

%% exact nonlinear operator W(x) (not regularised)
function z = nonLinearExact(x,curElem,lvl,p)

mu1 = p.problem.mu1;
mu2 = p.problem.mu2;
t1 = p.problem.t1;
t2 = p.problem.t2;

z = zeros(length(x),1);

index1 = find( x <= t1 );
z(index1) = mu2 * 1/2 * x(index1).^2;

index2 = find( x > t1 & x < t2);
z(index2) = t1*mu2* ( x(index2) - 1/2 * t1);

index3 = find( x >= t2 );
z(index3) = mu1 * 1/2 * x(index3).^2 + 1/2 * mu1 * t2^2 - 1/2 * mu2 * t1^2;


%% exact derivative of nonlinear operator, DW(x)/x (not regularised)
function z = nonLinearExactDer(x,curElem,lvl,p)

mu1 = p.problem.mu1;
mu2 = p.problem.mu2;
t1 = p.problem.t1;
t2 = p.problem.t2;

z = zeros(length(x),1);

index1 = find( x <= t1 );
z(index1) = mu2;

index2 = find( x > t1 & x < t2);
z(index2) = t1*mu2* 1/( x(index2));

index3 = find( x >= t2 );
z(index3) = mu1;

%% exact second derivative of nonlinear operator, D^2W(x) (not regularised)
function z = nonLinearExactSecDer(x,curElem,lvl,p)

mu1 = p.problem.mu1;
mu2 = p.problem.mu2;
t1 = p.problem.t1;
t2 = p.problem.t2;

z = zeros(length(x),1);

index1 = find( x <= t1 );
z(index1) = mu2;

index2 = find( x > t1 & x < t2);
z(index2) = 0;

index3 = find( x >= t2 );
z(index3) = mu1;
