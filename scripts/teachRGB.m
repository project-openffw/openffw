lineWidth = 2;

c4n = [-1 0;
		1 0;
		0 1;
		0 0;
		 1/2 1/2;
		 -1/2 1/2];

	
n4e = [1 2 3];

		 
% plot triangle
figure
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'LineWidth',lineWidth);
plot(c4n([1 2],1),c4n([1 2],2),'k','LineWidth',lineWidth+3);
axis off
hold off

% plot green
figure
hold on

trisurf([1 4 3;4 2 3],c4n(:,1),c4n(:,2),zeros(size(c4n,1),1),'FaceColor','g','LineWidth',lineWidth-.5);
plot(c4n([1 3],1),c4n([1 3],2),'k','LineWidth',lineWidth+2);
plot(c4n([3 2],1),c4n([3 2],2),'k','LineWidth',lineWidth+2);
axis off
hold off

% plot blueR
figure
hold on
trisurf([2 5 4;1 4 3;4 5 3],c4n(:,1),c4n(:,2),zeros(size(c4n,1),1),'FaceColor','b','LineWidth',lineWidth-.5);
plot(c4n([3 1],1),c4n([3 1],2),'k','LineWidth',lineWidth+2);
plot(c4n([3 4],1),c4n([3 4],2),'k','LineWidth',lineWidth+2);
plot(c4n([4 2],1),c4n([4 2],2),'k','LineWidth',lineWidth+2);
axis off
hold off

% plot red
figure
hold on
trisurf([1 4 6;4 5 6;4 2 5;5 3 6],c4n(:,1),c4n(:,2),zeros(size(c4n,1),1),'FaceColor','r','LineWidth',lineWidth-.5);
plot(c4n([1 4],1),c4n([1 4],2),'k','LineWidth',lineWidth+2);
plot(c4n([4 2],1),c4n([4 2],2),'k','LineWidth',lineWidth+2);
plot(c4n([5 6],1),c4n([5 6],2),'k','LineWidth',lineWidth+2);
axis off
hold off


