function p = OptimalDesign_Square_exact3(p)
%author: David Guenther
%Yosida regularization  
%given W*_epsilon by regularisation of W*
%given W_epsilon by conjugation of W*_epsilon
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

%% PDE definition
p.problem.geom = 'Square';

p.problem.epsilon = 1e-3;

p.problem.alpha = 0.5;

p = ODgenericNonLinear(p);
p = ODgetRegularConj(p);
% p = ODgetNonLinearRegularConj(p);

p.problem.u_D = @u_D;
p.problem.f = @f;
p.problem.sigma0 = @sigma0;
p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;

%% supply exact solution
function val = u_exact(points,p)

alpha = p.problem.alpha;

x = points(:,1);
y = points(:,2);

% val = zeros(length(x),1);
% 
% I1 = find(abs(x+y)>1.5);
% I2 = find(abs(x+y)<0.5);
% J = find(abs(x-y)>0.5);
% K = setdiff(1:length(x),[I1;I2;J]);
% 
% val(I1) = alpha*(abs(x(I1)+y(I1))-1);
% val(I2) = alpha*(abs(x(I2)+y(I2))-1);
% val(J) = alpha*(abs(x(J)-y(J))-1);
% val(K) = 0;
val = -alpha*(x+y);

%% supply gradient u_exact
function val = gradU_exact(points,p)

alpha = p.problem.alpha;

x = points(:,1);
y = points(:,2);

% val = zeros(length(x),2);

% I1 = find(abs(x+y)>1.5);
% I2 = find(abs(x+y)<0.5);
% J = find(abs(x-y)>0.5);
% K = setdiff(1:length(x),[I1;I2;J]);
% 
% val(I1,:) = alpha;
% val(I2,:) = alpha;
% val(J,1) = alpha; val(J,2) = -alpha;
% val(K,:) = 0;

val = -alpha*ones(length(x),2);

%% supply volume force f
function val = f(points,curElem,lvl,p)

val = zeros(length(points(:,1)),1);

%% supply Dirichlet boundary
function val = u_D(points,p)

% alpha = p.problem.alpha;
u_exact = p.problem.u_exact;
% val = alpha*x.*y;
val = u_exact(points,p);

%% supply exact stresses
function z = sigma0(points,curElem,lvl,p)

gradU_exact = p.problem.gradU_exact;
nonLinear = p.problem.nonLinearExactDer;

evalFunc = gradU_exact(points,p);
absFunc = (evalFunc(:,1).^2 + evalFunc(:,2).^2).^(1/2);
evalNonLinear = nonLinear(absFunc,curElem,lvl,p);

z = (evalNonLinear*[1,1]).*evalFunc;
