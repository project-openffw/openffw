function val = funcHandleStima(pts,curElem,curLvl,p)
% Function handle to calculate the local stiffness matrix.

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

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = p.problem.kappa;
gradBasis  = p.statics.gradBasis;

%% Calculate local stiffness matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

curKappa = kappa(pts,p);
curGrad = gradBasis(pts,curElem,curLvl,p);

nrBasisFuncU = size(curGrad,1);

nrPts = size(pts,1);
val = zeros(nrBasisFuncU,nrBasisFuncU,nrPts);

for curPt = 1:nrPts    
    val(:,:,curPt) = curGrad(:,:,curPt)*curKappa(:,:,curPt)'*curGrad(:,:,curPt)'; 
end