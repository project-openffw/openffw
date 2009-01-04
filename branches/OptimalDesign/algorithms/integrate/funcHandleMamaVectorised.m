function val = funcHandleMamaVectorised(pts,pts_ref,parts,curLvl,p)
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

%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = p.problem.mu;
basis = p.statics.basisVectorised;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrPts = size(pts,1);
curMu = mu(pts,p);
curBasis = basis(pts,pts_ref,parts,curLvl,p);
curBasis = reshape(curBasis',[size(curBasis,2) 1 size(curBasis,1)]);
curMu = reshape( (curMu(:)*ones(1,size(curBasis,1)^2))',[size(curBasis,1) size(curBasis,1) nrPts]);

val = matMul(curBasis,permute(curBasis,[ 2 1 3 ]));
val = curMu.*val;


