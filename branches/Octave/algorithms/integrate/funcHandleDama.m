function val = funcHandleDama(x,y,curElem,curLvl,p)
% Function handle to calculate the local dama matrix.

% Copyright 2007 Joscha Gedicke
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

%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = p.problem.lambda;
gradBasis  = p.statics.gradBasis;
basis = p.statics.basis;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curLambda = lambda(x,y,p);
curGrad = gradBasis(x,y,curElem,curLvl,p);
curBasis = basis(x,y,curElem,curLvl,p);

nrBasisFuncU = size(curGrad,1);

val = zeros(nrBasisFuncU,nrBasisFuncU,length(x));

for i = 1 : length(x)    
    val(:,:,i) = curGrad(:,:,i)*curLambda(i,:)'*curBasis(i,:);    
end

