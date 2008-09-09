function [x,w] = getGaussPoints(n)
% Gives n Gauss Point for the unite vector
% for Gauss-Legendere integration
%
% author: Joscha Gedicke

%Find n Gauss_Legendre Points for Intervall [-1,1]
gamma = [1 : n-1] ./ sqrt(4*[1 : n-1].^2 - ones(1,n-1) );
[V,D] = eig( diag(gamma,1) + diag(gamma,-1) );
x = diag(D);
w = 2*V(1,:).^2;

%linear map to Intervall [0,1]
x = .5 * x + .5;
w = .5 * w';