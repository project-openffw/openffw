function p = drawUminusUh(p,lvl)


if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end


clf;
pdeSolver = p.params.pdeSolver;
customDrawUminusUh = str2func([pdeSolver,'drawUminusUh']);

p = customDrawUminusUh(p,lvl);