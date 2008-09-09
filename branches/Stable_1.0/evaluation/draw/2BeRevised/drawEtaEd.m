function p = drawEtaEd(p,lvl)

lineWidth = 4;
nrColors = 256;

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
n4ed = p.level(lvl).enum.n4ed;
etaEd = p.level(lvl).etaEd;

markedEdges = p.level(lvl).markededges;
nrNodes = p.level(lvl).nrNodes;
nrEdges = p.level(lvl).nrEdges;

figure
map = colormap(jet(nrColors));
hold on
% colormap('default');
etaEd = etaEd./max(etaEd);
etaEd = etaEd *(nrColors);

triplot(n4e,c4n(:,1),c4n(:,2),'LineWidth',1,'Color','b');
for j = 1:nrEdges
	colorIdx = floor(etaEd(j));
	if(colorIdx > 0)
		plot3(c4n(n4ed(j,:),1),c4n(n4ed(j,:),2),[etaEd(j),etaEd(j)],'LineWidth',lineWidth,'Color',map(colorIdx,:));
	end
end
% colorbar
colorbar('location','SouthOutside','XTickLabel',...
    {'Low Error','','','High Error',''})
title(sprintf('Estimated Error on Edges (Level = %g)',lvl));
