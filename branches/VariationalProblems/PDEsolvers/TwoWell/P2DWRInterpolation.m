function val = P2DWRInterpolation(x,y,curElem,lvl,p)

% Copyright 2007 Joscha Gedicke, Lena Noack
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


%% P2 Interpolation


%% Output
%p.level(lvl).P2DWRu4e = @getI2Auh;
%
%function val = getI2Auh(x,y,lvl,curElem,p)

nrNodes = p.level(end).nrNodes;
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];

Auh = p.statics.Au_h;
evalAuh = Auh(coords(:,1),coords(:,2),curElem,lvl,p);

basisP2 = p.statics.basisP2;
evalBasisP2 = basisP2(x,y,curElem,lvl,p);

I2Au_h = (evalAuh * evalBasisP2')';
I2Auh = repmat(I2Au_h,length(x),1);

size(I2Auh)


dofU4e = [n4e,nrNodes+ed4e]; 
size(dofU4e)
val = I2Auh(dofU4e);