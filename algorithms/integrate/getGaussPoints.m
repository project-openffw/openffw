function [x,w] = getGaussPoints(n)
% Gives n Gauss points for the unite vector.
% Used for Gauss-Legendere integration.
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
gamma = [1 : n-1] ./ sqrt(4*[1 : n-1].^2 - ones(1,n-1) );
[V,D] = eig( diag(gamma,1) + diag(gamma,-1) );
x = diag(D);
w = 2*V(1,:).^2;

%% linear map to Intervall [0,1]
x = .5 * x + .5;
w = .5 * w';
