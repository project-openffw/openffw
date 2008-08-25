function p = P1P1init(p)
% makes available all necessary initial data,
% handels the discrete displacement,stress, ...

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

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.statics.basisU = @getDisplacementBasis;
p.statics.u_h = @getU_h;
p.statics.sigma_h = @getSigma_h;

% integration parameters 
%  -> up to which polynomial degree shall integration be exact?
p.params.integrationDegrees.createLinSys.Rhs = 9;

%% Basis Functions etc. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_h = getU_h(pts,curElem,lvl,p)
nrPts = size(pts,1);
u = p.level(lvl).u4e;
basisU = p.statics.basisU;
evalBasisU = basisU(pts,curElem,lvl,p);
nrBasisFuncU = size(evalBasisU,3);
curU = u(:,:,curElem);
curU = curU(:)';
curU = repmat(curU,2*nrPts,1);
curU = reshape(curU,[nrPts,2,nrBasisFuncU]);

u_h = evalBasisU.*curU;
u_hU = u_h(:,:,1:nrBasisFuncU/2);
u_hV = u_h(:,:,nrBasisFuncU/2+1:end);

u_hU = squeeze(sum(u_hU,3));
u_hV = squeeze(sum(u_hV,3));

u_h = [u_hU(:,1)'; u_hV(:,2)'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigma_h = getSigma_h(pts,curElem,lvl,p)
nrPts = size(pts,1);
sigma = p.level(lvl).sigma;
curSigma = sigma(:,:,curElem);
sigma_h = repmat(curSigma,[1,1,nrPts]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basisU = getDisplacementBasis(pts,curElem,lvl,p)
nrPts = size(pts,1);
x = pts(:,1);
y = pts(:,2);

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
area4e = p.level(lvl).enum.area4e;

nodes = n4e(curElem,:);
coords = c4n(nodes,:);
area = area4e(curElem);

P1 = coords(1,:);
P2 = coords(2,:);
P3 = coords(3,:);

x = x';
y = y';

b1 = 1/2/area*( (P2(2)-P3(2))*x + (P3(1)-P2(1))*y + P2(1)*P3(2)-P3(1)*P2(2) );
b2 = 1/2/area*( (P3(2)-P1(2))*x + (P1(1)-P3(1))*y + P3(1)*P1(2)-P1(1)*P3(2) );
b3 = 1/2/area*( (P1(2)-P2(2))*x + (P2(1)-P1(1))*y + P1(1)*P2(2)-P2(1)*P1(2) );

basisU = zeros(nrPts,2,6);
basisU(:,1,1) = b1;
basisU(:,1,2) = b2;
basisU(:,1,3) = b3;
basisU(:,2,4) = b1;
basisU(:,2,5) = b2;
basisU(:,2,6) = b3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
