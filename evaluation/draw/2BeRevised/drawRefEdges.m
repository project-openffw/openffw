function p = drawRefEdges(p,lvl)

lineWidth = 2;

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
n4ed = p.level(lvl).enum.n4ed;
nrNodes = p.level(lvl).nrNodes;
nrEdges = p.level(lvl).nrEdges;
nrElems = p.level(lvl).nrElems;

figure
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'LineWidth',lineWidth);
for j = 1:nrElems
	plot(c4n(n4e(j,[1,2]),1),c4n(n4e(j,[1,2]),2),'g','LineWidth',lineWidth);
end

title(sprintf('Reference Edges (Level = %g)',lvl));
