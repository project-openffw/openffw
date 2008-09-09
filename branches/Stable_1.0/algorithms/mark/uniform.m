function p = uniform(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uniform refinement
refineElemsBisec5 = false(nrElems,1);
refineEdges = true(nrEdges,1);

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).refineEdges = refineEdges';
p.level(end).refineElemsBisec5 = refineElemsBisec5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
