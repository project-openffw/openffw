
figure(1)
c4n = [0 0; 1 0; 1 1; 0 1; 0 0];
plot(c4n(:,1),c4n(:,2),'k-','Linewidth',8)

axis off

saveas(gcf,'geomSquare','epsc2');
saveas(gcf,'geomSquare','fig');

figure(2)

c4n = [0 0; 1 0; 1 1; -1 1; -1 -1; 0 -1; 0 0];
plot(c4n(:,1),c4n(:,2),'k-','Linewidth',8)

axis off

saveas(gcf,'geomLshape','epsc2');
saveas(gcf,'geomLshape','fig');

figure(3)

c4n = [0.92387953251129	0.38268343236509
0.38268343236509	0.92387953251129
-0.38268343236509	0.92387953251129
-0.92387953251129	0.38268343236509
-0.92387953251129	-0.38268343236509
-0.38268343236509	-0.92387953251129
0.38268343236509	-0.92387953251129
0.92387953251129	-0.38268343236509
0.92387953251129	0.38268343236509];

plot(c4n(:,1),c4n(:,2),'k-','Linewidth',8)

axis off

saveas(gcf,'geomOctagon','epsc2');
saveas(gcf,'geomOctagon','fig');

figure(4)

c4n = [0 0; 1 0; 1 1; -1 1; -1 -1; 1 -1; 1 0];
plot(c4n(:,1),c4n(:,2),'k-','Linewidth',8)

axis off

saveas(gcf,'geomSlit','epsc2');
saveas(gcf,'geomSlit','fig');
