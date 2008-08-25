function [x,w] = getConProdGaussPoints(n)
% Gauss Points for the Reference Triangle
%    T = {x+y| 0<=x<=1, 0<=y<=1, x+y<=1}
% Integrates polynomials up to degree 2n-1 exact
% using the Stroud Conical Product rule with n^2
% quadrature Points.
%
% input:  n - number of points
% output: x - Gauss points
%         w - Gauss weights

% Copyright 2007 Joscha Gedicke
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


%% Find n Gauss_Legendre Points for Intervall [-1,1]
gamma = (1 : n-1) ./ sqrt(4*(1 : n-1).^2 - ones(1,n-1) );
[V,D] = eig( diag(gamma,1) + diag(gamma,-1) );
r = diag(D);
a = 2*V(1,:).^2;% norm-factor -> int(1,-1,1)=2

%% Find n Gauss_Jacobi Points for Intervall [-1,1]
delta = -1./(4*(1 : n).^2-ones(1,n));
gamma = sqrt((2 : n).*(1 : n-1)) ./ (2*(2 : n)-ones(1,n-1));
[V,D] = eig( diag(delta)+diag(gamma,1)+diag(gamma,-1) );
s = diag(D);
b = 2*V(1,:).^2; % norm-factor -> int((1-x),-1,1)=2

%% linear map to Intervall [0,1]
% w(x) = 1 changes norm-factor from 2 to 1
r = .5 * r + .5;
a = .5 * a';
% w(x) = (1-x) changes norm-factor from 2 to 1/2
s = .5 * s + .5;
b = .25 * b';

%% conical product [ s_j , r_i(1-s_j) ]  a_i*b_j
s = repmat(s',n,1); s = s(:);
r = repmat(r,n,1);
x = [ s , r.*(ones(n^2,1)-s) ];
w = a*b';
w = w(:);
    
