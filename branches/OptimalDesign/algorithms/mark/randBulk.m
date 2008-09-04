function p = randBulk(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etaEd = loadField('p.level(end)','etaEd',p);

ed4e = p.level(end).enum.ed4e;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;
level = p.level(end).level;

thetaEd = loadField('p.params.modules.mark.bulk','thetaEd',p,1/2);

refineFirstLevel = loadField('p.params.modules.mark','refineFirstLevel',p,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bulk criterion
refineEdges = false(nrEdges,1);

if(level <= refineFirstLevel)
	refineEdges = true(nrEdges,1);
else	
	if( (nnz(etaEd) ~= 0) && level > 1)
        sumEtaEd = 0;
        while sumEtaEd < thetaEd * norm(etaEd,2)^2
            index = floor(nrEdges*rand(1));
            sumEtaEd = sumEtaEd + etaEd(index);
            etaEd(index) = 0;
            refineEdges(index) = true;
        end
    end
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).refineEdges = refineEdges';
p = closure(p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
