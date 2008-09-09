function val = funcHandleRHSVolumeVectorised(x,y,x_ref,y_ref,parts,curLvl,p)
% function handle for RHS f*phi_j
% use: val = funcHandleRHSVolume(x,y,x_ref,y_ref,parts,curlvl,p)

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
f = p.problem.f;
basisU = p.statics.basisVectorised;
dofU4e = p.level(curLvl).enum.dofU4e;
nrBasisFuncU = size(dofU4e,2);

%% f*phi_j %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
evalF = f(x,y,p);
dimf = size(evalF,2);
evalBasisU = basisU(x,y,x_ref,y_ref,parts,curLvl,p);
evalBasisU = reshape(evalBasisU, [length(x) dimf nrBasisFuncU]);
evalF = evalF(:)*ones(1,nrBasisFuncU);
evalF = reshape(evalF,[length(x) dimf nrBasisFuncU ]);
integrand = evalF .* evalBasisU;
integrand = sum(integrand,2);

val = permute(integrand, [ 3 2 1 ]);
