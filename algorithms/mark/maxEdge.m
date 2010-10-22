function p = maxEdge(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etaEd = loadField('p.level(end)','etaEd',p);

nrEdges = p.level(end).nrEdges;
level = p.level(end).level;

refineFirstLevel = loadField('p.params.modules.mark','refineFirstLevel',p,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximum criterion
refineEdges = false(nrEdges,1);

if(level <= refineFirstLevel)
	refineEdges = true(nrEdges,1);
else
	if( (nnz(etaEd) ~= 0) && level > 1)
        maxErrorEdge = find(max(etaEd) == etaEd);
		refineEdges(maxErrorEdge) = true;
    end
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).refineEdges = refineEdges';
p = closure(p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
