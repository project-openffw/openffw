function p = drawGrid(p,lvl)

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set graphic options from structure p
color = loadField('p.params.output','color',p,'k');
lineWidth = loadField('p.params.output','lineWidth',p,0.1);
drawInfo = loadField('p.params.output','drawInfo',p,false);
fontSize = loadField('p.params.output','fontSize',p,18);

% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

% load enumerated data
nrNodes = p.level(lvl).nrNodes;
nrElems = p.level(lvl).nrElems;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coordX = reshape(c4n(n4e',1),3,nrElems);
coordY = reshape(c4n(n4e',2),3,nrElems);
val = zeros(size(coordX));
patch(coordX,coordY,val,'FaceColor','none','LineWidth',lineWidth,'EdgeColor',color);

set(gca,'FontSize',fontSize);

if(drawInfo)
	xlabel(sprintf('Nr of Nodes %g',nrNodes));
	title(sprintf('Grid on Level %g',lvl));
end
