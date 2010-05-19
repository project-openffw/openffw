function p = ODP2getPhFuncs(p)
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

p.statics.sigma_h = @getSigma_h;

%% supply p_h = DW(\grad u_h)
function sigma_h = getSigma_h(points,curElem,lvl,p)

nonLinearExactDer = p.problem.nonLinearExactDer;
grad_h = p.statics.grad_h;

evalGrad = grad_h(points,curElem,lvl,p);
absGrad = ( evalGrad(:,1).^2 + evalGrad(:,2).^2 ).^(1/2);

sigma_h = (nonLinearExactDer(absGrad,curElem,lvl,p)*[1 1]).*evalGrad;
