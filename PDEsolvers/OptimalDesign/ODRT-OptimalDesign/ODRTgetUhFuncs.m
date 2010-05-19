function p = ODRTgetUhFuncs(p)
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
p.statics.Au_h = @getAu_h;
p.statics.gradAu_h = @getGradAu_h;
p.statics.I2u_h = @getI2u_h;
p.statics.gradI2u_h = @getGradI2u_h;

%% supply u_h
function u_h = getU_h(points,curElem,lvl,p)

u = p.level(lvl).u4e;
curU = u(curElem);
u_h = repmat(curU,length(points(:,1)),1);

%% supply Au_{h}
function Auh = getAu_h(points,curElem,lvl,p)

Auh4n = p.level(lvl).Auh;
n4e = p.level(lvl).geom.n4e;

basisP1 = getP1basis(points,curElem,lvl,p);

Auh = Auh4n(n4e(curElem,:),1)'*basisP1;

%% supply grad Au_{h}
function gradAuh = getGradAu_h(points,curElem,lvl,p)

Auh4n = p.level(lvl).Auh;
n4e = p.level(lvl).geom.n4e;

gradBasisP1 = getGradP1basis(points,curElem,lvl,p);

Auh = Auh4n(n4e(curElem,:),1)';

gradAuh = squeeze(matMul(repmat(Auh,[1 1 length(points(:,1))]),gradBasisP1))';

%% supply I^2(Au_h)
function I2Auh = getI2u_h(points,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];

Auh = p.statics.Au_h;
evalAuh = Auh(coords(:,1),coords(:,2),curElem,lvl,p);

basisP2 = p.statics.basisP2;
evalBasisP2 = basisP2(points,curElem,lvl,p);

I2Auh = (evalAuh * evalBasisP2)';

%% supply grad I^2(Au_h)
function gradI2uh = getGradI2u_h(points,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

Auh = p.statics.Au_h;
gradBasisP2 = p.statics.gradBasisP2;

coords = [c4n(n4e(curElem,:),:);
          midPoint4ed(ed4e(curElem,:),:)];
      
evalAuh = Auh(coords(:,1),coords(:,2),curElem,lvl,p);

evalGradBasisP2 = gradBasisP2(points,curElem,lvl,p);

evalAuh = repmat(evalAuh',[1 1 length(points(:,1))]);

gradI2uh = squeeze(matMul(evalAuh,evalGradBasisP2));
