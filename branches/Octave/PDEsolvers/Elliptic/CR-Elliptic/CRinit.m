function p = CRinit(p)
% makes available all necessary initial data,
% handels the discrete displacement,flux, ...

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

%% OUTPUT
p.statics.basisU = @getDisplacementBasis;
p.statics.u_h = @getU_h;
p.statics.sigma_h = @getSigma_h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_h = getU_h(x,y,curElem,lvl,p)

u = p.level(lvl).u4e;
basisU = p.statics.basisU;
basisU = basisU(x,y,curElem,lvl,p);

curU = u(curElem,:);
u_h = curU * basisU';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigma_h = getSigma_h(x,y,curElem,curKappa,lvl,p)

grad4e = p.level(lvl).grad4e;
curGrad = grad4e(curElem,:);

curGrad = reshape(curGrad'*ones(1,length(x)),[2 1 length(x)]);
sigma_h = matMul(curKappa,curGrad);
sigma_h = reshape(sigma_h,2,[])';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basisU = getDisplacementBasis(x,y,curElem,lvl,p)

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;
area4e = p.level(lvl).enum.area4e;

edges = ed4e(curElem,:);
coords = midPoint4ed(edges,:);
area = area4e(curElem);

area = area/4;

P1 = coords(1,:);
P2 = coords(2,:);
P3 = coords(3,:);

x = x';
y = y';

b1 = 1/2/area*( (P2(2)-P3(2))*x + (P3(1)-P2(1))*y + P2(1)*P3(2)-P3(1)*P2(2) );
b2 = 1/2/area*( (P3(2)-P1(2))*x + (P1(1)-P3(1))*y + P3(1)*P1(2)-P1(1)*P3(2) );
b3 = 1/2/area*( (P1(2)-P2(2))*x + (P2(1)-P1(1))*y + P1(1)*P2(2)-P2(1)*P1(2) );

basisU = [b1;b2;b3]';

