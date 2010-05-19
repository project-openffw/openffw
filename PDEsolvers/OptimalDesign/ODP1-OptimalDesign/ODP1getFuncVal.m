function p = ODP1getFuncVal(p)
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

%% INPUT
% load enumerated data
n4e = p.level(end).geom.n4e;
dofU4e = p.level(end).enum.dofU4e;
f4e = p.level(end).f4e;
lvl = size(p.level,2);
freeNodes = p.level(end).enum.freeNodes;
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,5);

%% compute the function value E(x)
ST = integrate(n4e,lvl,degree,@integrand,p);

I = dofU4e;
S = accumarray(I(:),ST(:));
rhs = accumarray(I(:),f4e(:));

funcVal = S - rhs;

%% OUTPUT 
p.level(end).funcVal = funcVal(freeNodes);

%% supply integrand: DW(\nabla u_h)*\nabla w_h
function val = integrand(points,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
stressBasis = p.statics.stressBasis;

x = points(:,1)';

evalSigma = sigma_h(points,curElem,lvl,p);
evalBasis = stressBasis(points,curElem,lvl,p);

evalSigma = reshape(evalSigma',[2 1 length(x)]);

val = matMul(evalBasis,evalSigma);
