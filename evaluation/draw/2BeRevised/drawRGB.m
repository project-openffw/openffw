function p = drawRGB(p,lvl)

lineWidth = 1;

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
ed4e = p.level(lvl).enum.ed4e;
n4ed = p.level(lvl).enum.n4ed;
newNode4ed = p.level(lvl).enum.newNode4ed;
markededges = p.level(lvl).markededges;

newNode4e = newNode4ed(ed4e);
unrefinedElems = find( all(newNode4e == 0 ,2) );
markedElems = find( any(newNode4e,2) );
nrMarkedEd4MarkedElems = sum(markededges( ed4e(markedElems,:) ),2);

I = find(nrMarkedEd4MarkedElems == 1);
gElems = markedElems(I);
I = find(nrMarkedEd4MarkedElems == 2);
bElems = markedElems(I);
I = find(nrMarkedEd4MarkedElems == 3);
rElems = markedElems(I);

figure
hold on

trisurf(n4e(gElems,:),c4n(:,1),c4n(:,2),zeros(size(c4n,1),1),'FaceColor','g');
trisurf(n4e(bElems,:),c4n(:,1),c4n(:,2),zeros(size(c4n,1),1),'FaceColor','b');
trisurf(n4e(rElems,:),c4n(:,1),c4n(:,2),zeros(size(c4n,1),1),'FaceColor','r');
triplot(n4e,c4n(:,1),c4n(:,2),'LineWidth',lineWidth,'Color','k');

for j = find(markededges)
	plot(c4n(n4ed(j,:),1),c4n(n4ed(j,:),2),'m','LineWidth',lineWidth+2);
end


hold off;

view(0,90)
title(sprintf('RGB (Level = %g)',lvl));
