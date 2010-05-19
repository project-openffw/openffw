function p = genericNonLinear(p)
% Copyright 2007 David Guenther
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%

p.problem.nonLinearConjExact = @nonLinearConjExact;
p.problem.nonLinearExact = @nonLinearExact;
p.problem.nonLinearExactDer = @nonLinearExactDer;
p.problem.nonLinearExactSecDer = @nonLinearExactSecDer;

%% exact conjugate nonlinear operator W*(x) (not regularised)
function z = nonLinearConjExact(x,curElem,lvl,p)

mu1 = p.problem.mu1;
mu2 = p.problem.mu2;
lambda = p.problem.lambda;
t1 = p.problem.t1;

z = zeros(length(x),1);

index1 = find( x <= mu2*t1 );
z(index1) = 1/2*x(index1).^2/mu2 + mu2*lambda;

index2 = find( x > mu2*t1 );
z(index2) = 1/2*x(index2).^2/mu1 + mu1*lambda;

%% exact nonlinear operator W(x) (not regularised)
function z = nonLinearExact(x,curElem,lvl,p)

mu1 = p.problem.mu1;
mu2 = p.problem.mu2;
lambda = p.problem.lambda;
t1 = p.problem.t1;
t2 = p.problem.t2;

z = zeros(length(x),1);

index1 = find( x <= t1 );
z(index1) = mu2*(1/2*x(index1).^2 - lambda);

index2 = find( x > t1 & x < t2);
z(index2) = t1*mu2*x(index2) - lambda*(mu1+mu2);

index3 = find( x >= t2 );
z(index3) = mu1*(1/2*x(index3).^2 - lambda);


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
z(index2) = t1*mu2./x(index2);

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
