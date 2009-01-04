function val = funcHandleNormalJumpVectorised(pts,pts_ref,parts,curLvl,p)
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
normals4ed = p.level(curLvl).enum.normals4ed(parts,:);
g          = p.problem.g;
NbEd       = p.level(curLvl).enum.NbEd;
DbEd       = p.level(curLvl).enum.DbEd;
gradU_h    = p.statics.gradU_hVectorised;
kappa      = p.problem.kappa;
midPoint4e = p.level(curLvl).enum.midPoint4e;

curElem   = e4ed(parts,:);

%% Jump of Fux In Normal Direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrPts = size(pts,1);
curMidPoint = midPoint4e(curElem(:,1),:);
curKappa = permute(kappa(curMidPoint,p),[3 2 1]);
curGradU_h = permute(gradU_h(pts,pts_ref,curElem(:,1),curLvl,p),[3 2 1]);

part_1 = (curKappa(:,1,1).*curGradU_h(:,1) + curKappa(:,1,2).*curGradU_h(:,2)).*normals4ed(:,1)+ ...
         (curKappa(:,2,1).*curGradU_h(:,1) + curKappa(:,2,2).*curGradU_h(:,2)).*normals4ed(:,2) ;

inner = find(curElem(:,2)>0);

curMidPoint = midPoint4e(curElem(inner,2),:);
curKappa = permute(kappa(curMidPoint,p),[3 2 1]);
curGradU_h = permute(gradU_h(pts(inner,:),pts_ref,curElem(inner,2),curLvl,p),[3 2 1]);

part_2 = zeros(size(part_1));
part_2(inner) = (curKappa(:,1,1).*curGradU_h(:,1) + curKappa(:,1,2).*curGradU_h(:,2)).*normals4ed(inner,1)+ ...
                (curKappa(:,2,1).*curGradU_h(:,1) + curKappa(:,2,2).*curGradU_h(:,2)).*normals4ed(inner,2) ;
     
     
[c, ia, nb] = intersect(NbEd, parts);
if ~isempty(nb)
    part_2(nb) = g(pts(nb,:),normals4ed(nb,:),p);
end
[c, ia, db] = intersect(DbEd, parts);
if ~isempty(db)
    part_2(db) = part_1(db);
end

%% Return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val(1,1,:) = (part_1 - part_2).^2;

