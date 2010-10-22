function p = drawU(p,lvl)


if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end


clf;
pdeSolver = p.params.pdeSolver;
customDrawU = str2func([pdeSolver,'drawU']);

p = customDrawU(p,lvl);