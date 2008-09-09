function p = drawCondNr(p,lvl,lineStyle,myColor,marker,minDoF,drawInfo)


if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

if(nargin < 3 || isempty(lineStyle))
	lineStyle = '-';
end


if(nargin < 4 || isempty(myColor))
	myColor = 'k';
end

if(nargin < 5 || isempty(marker))
	marker = 'x';
end

if(nargin < 6 || isempty(minDoF))
	minDoF = 1;
end

if(nargin < 7 || isempty(drawInfo))
	drawInfo = false;
end

nrLevels = p.level(end).level;
refineFirstLevel = p.params.modules.mark.refineFirstLevel;
p = conditionNumber(p);

for curLvl = refineFirstLevel+1:lvl
	nrDoF4lvl(curLvl-1) = p.level(curLvl).nrDoF;
	cond4lvl(curLvl-1) = p.level(curLvl).conditionNr;
end

% lambda = (p.PDE.lambda)^(1/3);

I = find(nrDoF4lvl > minDoF);
nrDoF4lvl = nrDoF4lvl(I);
cond4lvl = cond4lvl(I);
% cond4lvl = cond4lvl(I)/lambda;

if(isempty(nrDoF4lvl))
	return
end

	loglog(nrDoF4lvl(1:end),cond4lvl,'Marker',marker,'linestyle',lineStyle,'Color',myColor);
	set(gca,'XScale','log');
	set(gca,'YScale','log');
	
	if(drawInfo)
		xlabel(sprintf('Nr of degrees of freedom'));
		title('Condition number of global energy matrix')
	end
	grid on
end
