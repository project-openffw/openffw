function p = maximum(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etaT = loadField('p.level(end)','etaT',p);
etaEd = loadField('p.level(end)','etaEd',p);
etaOsc = loadField('p.level(end)','etaOsc',p);

ed4e = p.level(end).enum.ed4e;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;
level = p.level(end).level;

thetaEd = loadField('p.params.modules.mark.bulk','thetaEd',p,1/2);
thetaT = loadField('p.params.modules.mark.bulk','thetaT',p,1/2);
thetaOsc = loadField('p.params.modules.mark.bulk','thetaOsc',p,1/2);

refineFirstLevel = loadField('p.params.modules.mark','refineFirstLevel',p,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maximum criterion
refineEdges = false(nrEdges,1);
refineElems = false(nrElems,1);
refineElemsBisec5 = false(nrElems,1);

if(level <= refineFirstLevel)
	refineEdges = true(nrEdges,1);
else
	if( (nnz(etaEd) ~= 0) && level > 1)
		refineEdges = (etaEd > thetaEd * max(etaEd));
	end

	if( (nnz(etaT) ~= 0) && level > 1)
        I = (etaT > thetaT * max(etaT))';
		refineElems(I) = true;
    end
    
    if( (nnz(etaOsc) ~= 0) && level > 1)
        I = (etaOsc > thetaOsc * max(etaOsc))';
		refineElemsBisec5(I) = true;	
    end
    
   refineEdges4e = ed4e((refineElems | refineElemsBisec5),:);
   refineEdges(refineEdges4e(:)) = true;
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).refineEdges = refineEdges';
p.level(end).refineElemsBisec5 = refineElemsBisec5';
p = closure(p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
