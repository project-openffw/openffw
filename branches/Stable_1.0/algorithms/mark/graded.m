function p = graded(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
level = p.level(end).level;

beta = loadField('p.params.modules.mark.graded','beta',p,1/3);
N = loadField('p.params.modules.mark.graded','gradeN',p,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = (1/2)^N;
% max_h = (6/7)^N;
mu = beta;


count = 1;
nextLevel = level + 1;

p.level(nextLevel).geom = p.level(nextLevel-1).geom;
enumerate = p.statics.enumerate;
refine = p.statics.refine;
 p = enumerate(p);
% 	p.level(nextLevel-1).enum.newNode4ed = [];

while count ~= 0
    count = 0;
    
    ed4e = p.level(end).enum.ed4e;
    midPoint4e = p.level(end).enum.midPoint4e;
    area4e = p.level(end).enum.area4e;
    nrEdges = p.level(end).nrEdges;
    nrElems = p.level(end).nrElems;
    
    
    refineEdges = false(nrEdges,1);
    refineElems = false(nrElems,1);
    
    for curElem = 1:nrElems
        mP = midPoint4e(curElem,:);
        area = area4e(curElem,:);
        hT = sqrt(area);
        PhiMuT = norm(([0,0] - mP))^(1-mu);
        if hT > H*PhiMuT;
            refineElems(curElem) = true;
            count = count + 1;
        end
    end

    refineEdges4e = ed4e(refineElems,:);
    refineEdges(refineEdges4e(:)) = true;
    
    p.level(end).refineEdges = refineEdges';
    q = closure(p);
    q.level(end).refineElemsBisec5 = false(nrElems,1);
    q = refine(q);
    p.level(end).geom = q.level(end).geom;
    p = enumerate(p);
end 

nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;
refineEdges = false(nrEdges,1);

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).refineEdges = refineEdges';
p.level(end).refineElemsBisec5 = false(nrElems,1);
p = closure(p);
p.params.modules.mark.graded.gradeN = N+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%