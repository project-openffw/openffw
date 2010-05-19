function p = OptimalDesign_Octagon(p)
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
p.problem.geom = 'Octagon';

p.problem.epsilon = 1e-3;

p = ODgenericNonLinear(p);
p = ODgetRegularConj(p);
p = ODgetNonLinearRegularConj(p);

p.problem.u_D = @u_D;
p.problem.f = @f;

%% Dirichlet data
function val = u_D(points,p)

val = zeros(length(points(:,1)),1);

%% Volume force
function val = f(points,curElem,lvl,p)

val = ones(length(points(:,1)),1);
