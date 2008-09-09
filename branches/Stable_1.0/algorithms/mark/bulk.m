function p = bulk(p)

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

% Bulk criterion
refineEdges = false(nrEdges,1);
refineElems = false(nrElems,1);
refineElemsBisec5 = false(nrElems,1);

if(level <= refineFirstLevel)
	refineEdges = true(nrEdges,1);
else	
	if( (nnz(etaEd) ~= 0) && level > 1)
		[sortedEtaEd,I] = sort(etaEd,'descend');
		sumEtaEd = cumsum(sortedEtaEd.^2);
		k = find(sumEtaEd >= thetaEd * norm(etaEd,2)^2,1,'first');
		
		if(thetaEd == 1)
			refineEdges = true(nrEdges,1);
		else
			refineEdges(I(1:k)) = true;
		end
	end
	
	if( (nnz(etaT) ~= 0) && level > 1)
		[sortedEtaT,I] = sort(etaT,'descend');
		sumEtaT = cumsum(sortedEtaT.^2);
		k = find(sumEtaT >= thetaT * norm(etaT,2)^2,1,'first');
		refineElems(I(1:k)) = true;
	end
	
	if( (nnz(etaOsc) ~= 0) && level > 1)
		[sortedEtaOsc,I] = sort(etaOsc,'descend');
		sumEtaOsc = cumsum(sortedEtaOsc.^2);
		k = find(sumEtaOsc >= thetaOsc * norm(etaOsc,2)^2,1,'first');
		refineElemsBisec5(I(1:k)) = true;
	end
	
	refineEdges4e = ed4e((refineElems | refineElemsBisec5),:);
	refineEdges(refineEdges4e(:)) = true;
end

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).refineEdges = refineEdges';
p.level(end).refineElemsBisec5 = refineElemsBisec5;
p = closure(p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
