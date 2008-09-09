function p = STUBestimate(p)

nrEdges = p.level(end).nrEdges;

p.level(end).etaEd = zeros(nrEdges,1);
p.level(end).estimatedError = 0;
