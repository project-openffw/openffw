function val = funcHandleDamaVectorised(pts,pts_ref,parts,curLvl,p)
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

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = p.problem.lambda;
gradBasis  = p.statics.gradBasisVectorised;
basis = p.statics.basisVectorised;

%% Calculate local dama %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curLambda = lambda(pts,p);
curGrad = gradBasis(pts,pts_ref,parts,curLvl,p);
curBasis = basis(pts,pts_ref,parts,curLvl,p);

nrPts = size(pts,1);
dim = size(pts,2);
curLambda = reshape( curLambda',[dim 1 nrPts] );
curBasis = reshape(curBasis',[1 size(curBasis,2) size(curBasis,1)]);

val = matMul(curGrad,curLambda);
val = matMul(val,curBasis);