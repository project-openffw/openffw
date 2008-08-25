function p = Elliptic_Lshape_exact(p)
% example on the L shaped domain
%
%   -laplace(u) = 0 in Omega
%            u  = 0 on Db = conv{(0,0),(1,0)} U conv{(0,0),(0,-1)}
%        du/dn  = g(x,y) on Nb = (boundary of Gamma)\Db
%
% with known solution (in poolar coordiantes):
%     u = r^(2/3)sin(2/3*phi)

% Copyright 2007 Jan Reininghaus, David Guenther
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


% PDE definition
p.problem.geom = 'LshapeNeumann';%'LshapeArbitrary';%'LshapeHangingNodes';%
p.problem.f = @f;
p.problem.g = @g;
p.problem.u_D = @u_D;
p.problem.kappa = @kappa;
p.problem.Dkappa = @Dkappa;
p.problem.lambda = @lambda;
p.problem.mu = @mu;

% load exact solution
p.problem.u_exact = @u_exact;
p.problem.gradU_exact = @gradU_exact;
return


%% Volume force %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = f(pts,p)
nrPts = size(pts,1);
z = zeros(nrPts,1);

%% Dirichlet boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = u_D(pts,p)
z = u_exact(pts,p);

%% Neumann boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = g(pts,normals,p)
nrPts = size(pts,1);
[phi,r] = cart2pol(pts(:,1),pts(:,2));
phi = 2*pi*(phi<-eps)+phi; % This yields 0 <= phi < 2 pi
z = zeros(nrPts,1);
for j = 1:nrPts
    rotNormal = [cos(phi(j)) sin(phi(j));-sin(phi(j)) cos(phi(j))]*normals(j,:)';
    z(j) = 2/3*r(j)^(-1/3)*[sin(2/3*phi(j)), cos(2/3*phi(j))]*rotNormal;
end

%% elliptic PDE coefficent kappa ( div(kappa*grad_u) ) %%%%%%%%%%%%%%%%%%%%
function z = kappa(pts,p)
nrPts = size(pts,1);
dim = size(pts,2);
z = zeros(dim,dim,nrPts);
for curPt = 1:nrPts 
    z(:,:,curPt) = eye(dim);
end

function z = Dkappa(pts,p)
nrElems = size(pts,1);
% 3rd dimension is derivative with respect to space-dimensions
z = zeros(2,2,2,nrElems);

%% elliptic PDE coefficent lambda ( lambda*grad_u ) %%%%%%%%%%%%%%%%%%%%%%%
function z = lambda(pts,p)
nrPts = size(pts,1);
dim = size(pts,2);
z = zeros(nrPts,dim);
 
%% elliptic PDE coefficent mu ( mu*u ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = mu(pts,p)
nrPts = size(pts,1);
z = zeros(nrPts,1);
 
%% exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = u_exact(pts,p)
[phi,r] = cart2pol(pts(:,1),pts(:,2));
phi( phi<0 ) = phi( phi<0 ) + 2*pi;
z =  r.^(2/3).*sin(2/3*phi);

%% gradient of exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = gradU_exact(pts,p)
nrPts = size(pts,1);
dim = size(pts,2);
[phi,r] = cart2pol(pts(:,1),pts(:,2));
phi = 2*pi*(phi<-eps)+phi; % This yields 0 <= phi < 2 pi
z = zeros(nrPts,dim);
for j = 1:nrPts
    rotNormal = [cos(phi(j)) sin(phi(j));-sin(phi(j)) cos(phi(j))];
    z(j,:) = 2/3*r(j)^(-1/3)*[sin(2/3*phi(j)), cos(2/3*phi(j))]*rotNormal;
end
