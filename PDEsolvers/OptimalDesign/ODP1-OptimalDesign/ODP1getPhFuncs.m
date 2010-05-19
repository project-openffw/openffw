function p = ODP1getPhFuncs(p)
%author: David Guenther
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

p.statics.sigma_h = @getSigma_h;
p.statics.Ap_h = @getAp_h;

%% supply p_h = DW(\grad u_h)
function sigma_h = getSigma_h(points,curElem,lvl,p)

nonLinearExactDer = p.problem.nonLinearExactDer;

grad4e = p.level(lvl).grad4e;
curGrad = grad4e(curElem,:);

sigma_h = nonLinearExactDer(norm(curGrad),curElem,lvl,p)*curGrad;
sigma_h = ones(length(points(:,1)),1) * sigma_h;

%% supply average p_h: A(p_h)
function Aph = getAp_h(points,curElem,lvl,p)

Aph4n = p.level(lvl).Aph;
n4e = p.level(lvl).geom.n4e;
basisU = p.statics.basisU;
basisP1 = basisU(points,curElem,lvl,p)';

AphX = Aph4n(n4e(curElem,:),1)'*basisP1;
AphY = Aph4n(n4e(curElem,:),2)'*basisP1;

Aph = [AphX',AphY'];

