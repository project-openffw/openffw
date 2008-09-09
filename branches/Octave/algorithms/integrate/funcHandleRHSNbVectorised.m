function val = funcHandleRHSNbVectorised(x,y,x_ref,y_ref,parts,curLvl,p)
% function handle for g*phi_j
% use: val = funcHandleRHSNb(x,y,x_ref,y_ref,parts,curLvl,p)

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
g = p.problem.g;
basis = p.statics.basisVectorised;
dofU4e = p.level(curLvl).enum.dofU4e;
nrBasisFunc = size(dofU4e,2);
e4ed = p.level(curLvl).enum.e4ed;
curNormals4ed = p.level(curLvl).enum.normals4ed(parts,:);


%% g*phi_j%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
evalG = g(x,y,curNormals4ed,p);
dimg = size(evalG,2);
evalBasis = basis(x,y,x_ref,y_ref,e4ed(parts),curLvl,p);
evalBasis = reshape(evalBasis, [length(x) dimg nrBasisFunc]);
evalG = evalG(:)*ones(1,nrBasisFunc);
evalG = reshape(evalG,[length(x) dimg nrBasisFunc ]);
integrand = evalG .* evalBasis;
integrand = sum(integrand,2);

%% Return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val = permute(integrand, [ 3 2 1 ]);
