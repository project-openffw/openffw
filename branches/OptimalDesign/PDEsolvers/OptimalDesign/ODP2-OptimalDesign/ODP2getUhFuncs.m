function p = ODP2getUhFuncs(p)
% author: David Guenther 
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

p.statics.u_h = @getU_h;
p.statics.grad_h = @getGrad_h;

%% supply u_h
function u_h = getU_h(points,curElem,lvl,p)

u = p.level(lvl).u4e;
basisU = p.statics.basisU;
basisU = basisU(points,curElem,lvl,p);

curU = u(curElem,:);
u_h = (curU * basisU')';

%% supply \grad u_h
function grad_h = getGrad_h(points,curElem,lvl,p)

stressBasis = p.statics.stressBasis;
u = p.level(lvl).u4e;

evalBasis = stressBasis(points,curElem,lvl,p);
u = ones(length(points(:,1)),1)*u(curElem,:);
u = reshape(u',[1 6 length(points(:,1))]);

grad_h = squeeze(matMul(u,evalBasis))';
