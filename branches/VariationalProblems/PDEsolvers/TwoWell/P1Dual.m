function p = P1Dual(p)
% Copyright 2007 Joscha Gedicke, Lena Noack

nrEdges = p.level(end).nrEdges;

p.level(end).etaEd = zeros(nrEdges,1);
p.level(end).estimatedError = 0;