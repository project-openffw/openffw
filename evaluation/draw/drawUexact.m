function p = drawUexact(p,lvl)


if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end


clf;
pdeSolver = p.params.pdeSolver;
customDrawUexact = str2func([pdeSolver,'drawUexact']);

p = customDrawUexact(p,lvl);