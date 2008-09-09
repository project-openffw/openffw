function p = drawMarkedGrid(p,lvl)

lineWidth = 2;


n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
n4ed = p.level(lvl).enum.n4ed;
markedEdges = p.level(lvl).markededges;
nrNodes = p.level(lvl).nrNodes;
nrEdges = p.level(lvl).nrEdges;

figure
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'LineWidth',lineWidth);
for j = find(markedEdges)
	plot(c4n(n4ed(j,:),1),c4n(n4ed(j,:),2),'r','LineWidth',lineWidth);
end

title(sprintf('Marked Edges after Closure (Level = %g)',lvl));
