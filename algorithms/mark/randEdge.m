function p = randEdge(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrEdges = p.level(end).nrEdges;
level = p.level(end).level;

refineFirstLevel = loadField('p.params.modules.mark','refineFirstLevel',p,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximum criterion
refineEdges = false(nrEdges,1);

if(level <= refineFirstLevel)
	refineEdges = true(nrEdges,1);
else
	if level > 1 
        index = floor(nrEdges*rand(1));
		refineEdges(index) = true;
    end
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).refineEdges = refineEdges';
p = closure(p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
