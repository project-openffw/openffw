function p = P1Averaging(p,curLvl)

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

if nargin < 2
    curLvl = length(p.level);
end

%% Input
dofU4e = p.level(curLvl).enum.dofU4e;
patch4dof = p.level(curLvl).enum.area4n;
n4e = p.level(curLvl).geom.n4e;
degree = 1;

%% P1Averaging
%fluxRecovery = integrateVectorised(n4e,curLvl,degree,p.statics.sigma_hVectorised,p);
fluxRecovery = integrate(n4e,curLvl,degree,p.statics.sigma_h,p);
dummy1 = fluxRecovery(:,1)*ones(1,size(dofU4e,2));
dummy2 = fluxRecovery(:,2)*ones(1,size(dofU4e,2));
dummy1 = accumarray(dofU4e(:),dummy1(:))./patch4dof;
dummy2 = accumarray(dofU4e(:),dummy2(:))./patch4dof;
fluxRecovery = [dummy1,dummy2];


fluxRecovery4e(:,:,1) = dummy1(dofU4e);
fluxRecovery4e(:,:,2) = dummy2(dofU4e);
%size 24x3x2 ok

%% Output
p.level(curLvl).fluxRecovery = fluxRecovery;
p.level(curLvl).fluxRecovery4e = fluxRecovery4e;

p.statics.fluxRecoveryUh = @getFluxRecoveryUh;
p.statics.DivfluxRecoveryUh = @getDivfluxRecoveryUh;

%% Function Handles
function val = getFluxRecoveryUh(pts,pts_ref,parts,lvl,p)

fr4e1 = p.level(lvl).fluxRecovery4e(parts,:,1);
fr4e2 = p.level(lvl).fluxRecovery4e(parts,:,2); 
basisU = p.statics.basisU;

basisU = basisU(pts,pts_ref,parts,lvl,p);

val = [fr4e1 * basisU',fr4e2 * basisU'];
val = reshape(val,[2 1 size(basisU,1)]);

function val = getDivfluxRecoveryUh(pts,pts_ref,parts,lvl,p)

fr4e1 = p.level(lvl).fluxRecovery4e(parts,:,1);
fr4e2 = p.level(lvl).fluxRecovery4e(parts,:,2); 
basisGradU = p.statics.stressBasis;

basisGradU = basisGradU(pts,pts_ref,parts,lvl,p);
basisGradU = permute(basisGradU,[ 3 1 2]);

val = sum(fr4e1 .* basisGradU(:,:,1),2)+sum(fr4e2 .* basisGradU(:,:,2),2);
