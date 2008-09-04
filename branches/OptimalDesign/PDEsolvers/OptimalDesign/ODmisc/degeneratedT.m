
c4n = [0 0; 1 0; 0.5 1];
n4e = [1 2 3];

figure(1)
clf
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'Color','k','LineWidth',3)
text(0.5,-0.05,'E','Fontsize',40)
axis off

figure(2)
clf
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'Color','k','LineWidth',3)
triplot(n4e,[0 0.5 0.5],[0 0 1],'Color','k','LineWidth',3)
text(0.25,-0.05,'E','Fontsize',40)
axis off


figure(3)
clf
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'Color','k','LineWidth',3)
triplot(n4e,[0 0.5 0.5],[0 0 1],'Color','k','LineWidth',3)
triplot(n4e,[0 0.25 0.5],[0 0 1],'Color','k','LineWidth',3)
text(0.125,-0.05,'E','Fontsize',40)
axis off


figure(4)
clf
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'Color','k','LineWidth',3)
triplot(n4e,[0 0.5 0.5],[0 0 1],'Color','k','LineWidth',3)
triplot(n4e,[0 0.25 0.5],[0 0 1],'Color','k','LineWidth',3)
triplot(n4e,[0 0.125 0.5],[0 0 1],'Color','k','LineWidth',3)
text(0.0625,-0.05,'E','Fontsize',40)
axis off
