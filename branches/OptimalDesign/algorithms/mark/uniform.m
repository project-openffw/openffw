function p = uniform(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;
refineFirstLevel = loadField('p.params.modules.mark','refineFirstLevel',p,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uniform refinement
refineElemsBisec5 = false(nrElems,1);
if refineFirstLevel == 0
    refineEdges = false(nrEdges,1);
else
    refineEdges = true(nrEdges,1);
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).refineEdges = refineEdges';
p.level(end).refineElemsBisec5 = refineElemsBisec5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
