function val = funcHandleNormalJump(x,y,curEd,curLvl,p)

% author: Joscha Gedicke

%% INPUT
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

%% JUMP OF FLUX IN NORMAL DIRECTION
curMidPoint = midPoint4e(curElem(1),:);
curKappa = kappa(curMidPoint(1),curMidPoint(2),p);
part_1 = (gradU_h(x,y,curElem(1),curLvl,p)*curKappa')*curNormal';
if( curElem(2) > 0 )
    curMidPoint = midPoint4e(curElem(1),:);
    curKappa = kappa(curMidPoint(1),curMidPoint(2),p);
    part_2 = (gradU_h(x,y,curElem(2),curLvl,p)*curKappa')*curNormal';
elseif( isempty(Nb) || isempty(find(NbEd==curEd)) )
    part_2 = part_1;
else
    part_2 = g(x,y,(curNormal'*ones(1,length(x)))',p);
end

%% RETURN
val(1,:,:) = ((part_1 - part_2).^2)';
