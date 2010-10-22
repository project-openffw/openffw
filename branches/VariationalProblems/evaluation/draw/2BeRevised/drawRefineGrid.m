function p = drawRefineGrid(p,lvl)

lineWidth = 2;

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
n4ed = p.level(lvl).enum.n4ed;
refineEdges = p.level(lvl).refineEdges;
nrNodes = p.level(lvl).nrNodes;
nrEdges = p.level(lvl).nrEdges;

figure
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'LineWidth',lineWidth);
for j = find(refineEdges)
	plot(c4n(n4ed(j,:),1),c4n(n4ed(j,:),2),'r','LineWidth',lineWidth);
end

title(sprintf('Marked Edges before Closure (Level = %g)',lvl));
