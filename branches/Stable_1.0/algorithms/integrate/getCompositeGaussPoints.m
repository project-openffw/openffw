function [x,w] = getCompositeGaussPoints(n)
% Gives n Gauss Point for the unite square
% for Gauss-Legendere integration
%
% author: Joscha Gedicke

%Find n Gauss_Legendre Points for Intervall [-1,1]
gamma = [1 : n-1] ./ sqrt(4*[1 : n-1].^2 - ones(1,n-1) );
[V,D] = eig( diag(gamma,1) + diag(gamma,-1) );
x = diag(D);
w = 2*V(1,:).^2;

%linear map to square [0,1] times [0,1]
x1 = .5*x + .5;
x2 = .5*x + .5;
xi1 = repmat(x1',n,1);
x = [xi1(:),repmat(x2,n,1)];
w = w'*w;
w = .25 *w(:);
