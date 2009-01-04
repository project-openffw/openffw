function val = funcHandleMama(pts,curElem,curLvl,p)
% Function handle to calculate the local mass matrix.

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
mu = p.problem.mu;
basis = p.statics.basis;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curMu = mu(pts,p);
curBasis = basis(pts,curElem,curLvl,p);

nrBasisFuncU = size(curBasis,2);
nrPts = size(pts,1);
val = zeros(nrBasisFuncU,nrBasisFuncU,nrPts);

for curPt = 1:nrPts   
    val(:,:,curPt) = curMu(curPt)*curBasis(curPt,:)'*curBasis(curPt,:);    
end

