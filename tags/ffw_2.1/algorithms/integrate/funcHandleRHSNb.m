function val = funcHandleRHSNb(pts,curEd,curLvl,p)
% function handle for g*phi_j
% use: val = funcHandleRHSNb(x,y,curElem,curlvl,p)

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
g = p.problem.g;
basisU = p.statics.basisU;
dofU4e = p.level(curLvl).enum.dofU4e;
nrBasisFuncU = size(dofU4e,2);
e4ed = p.level(curLvl).enum.e4ed;
curNormal = p.level(curLvl).enum.normals4ed(curEd,:);


%% g*phi_j %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrPts = size(pts,1);
evalG = g(pts,(curNormal(:)*ones(1,nrPts))',p);
dimg = size(evalG,2);
evalBasisU = basisU(pts,e4ed(curEd),curLvl,p);
evalBasisU = reshape(evalBasisU, [nrPts dimg nrBasisFuncU]);
evalG = evalG(:)*ones(1,nrBasisFuncU);
evalG = reshape(evalG,[nrPts dimg nrBasisFuncU ]);
integrand = evalG .* evalBasisU;
integrand = sum(integrand,2);

%% Return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val = permute(integrand, [ 3 2 1 ]);
