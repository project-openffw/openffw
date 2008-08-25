function val = funcHandleNormalJump(pts,curEd,curLvl,p)
% function handle to calculate the jump in normal direction

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
e4ed       = p.level(curLvl).enum.e4ed;
normals4ed = p.level(curLvl).enum.normals4ed;
g          = p.problem.g;
Nb         = p.level(curLvl).geom.Nb;
NbEd       = p.level(curLvl).enum.NbEd;
gradU_h    = p.statics.gradU_h;
kappa      = p.problem.kappa;
midPoint4e = p.level(curLvl).enum.midPoint4e;

curElem   = e4ed(curEd,:);
curNormal = normals4ed(curEd,:);

%% Jump of Flux In Normal Direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrPts = size(pts,1);
curMidPoint = midPoint4e(curElem(1),:);
curKappa = kappa(curMidPoint(1),curMidPoint(2),p);
curGradU_h = permute(gradU_h(pts,curElem(1),curLvl,p),[3 2 1]);
part_1 = (curGradU_h*curKappa')*curNormal';
if( curElem(2) > 0 )
    curMidPoint = midPoint4e(curElem(2),:);
    curKappa = kappa(curMidPoint,p);
    curGradU_h = permute(gradU_h(pts,curElem(2),curLvl,p),[3 2 1]);
    part_2 = (curGradU_h*curKappa')*curNormal';
elseif( isempty(Nb) || isempty(find(NbEd==curEd, 1)) )
    part_2 = part_1;
else
    part_2 = g(pts,(curNormal'*ones(1,nrPts))',p);
end

%% Return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val(1,:,:) = ((part_1 - part_2).^2)';
