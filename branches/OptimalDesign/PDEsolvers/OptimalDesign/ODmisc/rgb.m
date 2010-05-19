% Copyright 2007 David Guenther
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%


c4n = [0 0; 1 0; 0 1];
n4e = [1 2 3];

figure(1)
clf
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'Color','k','LineWidth',2)
line([1 0],[0 1],'Linewidth',8,'Color','k')
text(0.55,0.55,'E_{ref}','Fontsize',30)
axis off
saveas(gcf,'initalT','epsc2');
saveas(gcf,'initalT','fig');

figure(2)
clf
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'Color','k','LineWidth',2)
line([1 0],[0 1],'Linewidth',3,'Color','k')
y = 1 - [0:0.1:1];
plot([0:0.1:1],y,'k*','Markersize',20)
axis off
saveas(gcf,'greenMarkT','epsc2');
saveas(gcf,'greenMarkT','fig');

figure(3)
clf
hold on
n4eNew = [1 2 4; 1 4 3];
c4nNew = [c4n;[0.5 0.5]]
triplot(n4eNew,c4nNew(:,1),c4nNew(:,2),'Color','k','LineWidth',2)
line([0 1],[0 0],'Linewidth',8,'Color','k')
line([0 0],[0 1],'Linewidth',8,'Color','k')
axis off
saveas(gcf,'greenT','epsc2');
saveas(gcf,'greenT','fig');

figure(4)
clf
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'Color','k','LineWidth',2)
line([1 0],[0 1],'Linewidth',3,'Color','k')
y = 1 - [0:0.1:1];
plot([0:0.1:1],y,'k*','Markersize',20)
plot([0:0.1:1],0*y,'k*','Markersize',20)
axis off
saveas(gcf,'blueLMarkT','epsc2');
saveas(gcf,'blueLMarkT','fig');

figure(5)
clf
hold on
n4eNew = [1 5 4; 5 2 4;1 4 3];
c4nNew = [c4n;[0.5 0.5;0.5 0]]
triplot(n4eNew,c4nNew(:,1),c4nNew(:,2),'Color','k','LineWidth',2)
line([0 0],[0 1],'Linewidth',8,'Color','k')
line([0 0.5],[0 0.5],'Linewidth',8,'Color','k')
line([0.5 1],[0.5 0],'Linewidth',8,'Color','k')
axis off
saveas(gcf,'blueLT','epsc2');
saveas(gcf,'blueLT','fig');

figure(6)
clf
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'Color','k','LineWidth',2)
line([1 0],[0 1],'Linewidth',3,'Color','k')
y = 1 - [0:0.1:1];
plot([0:0.1:1],y,'k*','Markersize',20)
plot(0*y,[0:0.1:1],'k*','Markersize',20)
axis off
saveas(gcf,'blueRMarkT','epsc2');
saveas(gcf,'blueRMarkT','fig');

figure(7)
clf
hold on
n4eNew = [1 2 4; 1 4 5;5 4 3];
c4nNew = [c4n;[0.5 0.5;0 0.5]];
triplot(n4eNew,c4nNew(:,1),c4nNew(:,2),'Color','k','LineWidth',2)
line([0 1],[0 0],'Linewidth',8,'Color','k')
line([0 0.5],[0 0.5],'Linewidth',8,'Color','k')
line([0.5 0],[0.5 1],'Linewidth',8,'Color','k')
axis off
saveas(gcf,'blueRT','epsc2');
saveas(gcf,'blueRT','fig');

figure(8)
clf
hold on
triplot(n4e,c4n(:,1),c4n(:,2),'Color','k','LineWidth',2)
line([1 0],[0 1],'Linewidth',3,'Color','k')
y = 1 - [0:0.1:1];
plot([0:0.1:1],y,'k*','Markersize',20)
plot([0:0.1:1],0*y,'k*','Markersize',20)
plot(0*y,[0:0.1:1],'k*','Markersize',20)
axis off
saveas(gcf,'redMarkT','epsc2');
saveas(gcf,'redMarkT','fig');

figure(9)
clf
hold on
n4eNew = [1 4 6; 4 2 5;6 5 3;4 5 6];
c4nNew = [c4n;[0.5 0;0.5 0.5;0 0.5]];
triplot(n4eNew,c4nNew(:,1),c4nNew(:,2),'Color','k','LineWidth',2)
line([0 1],[1 0],'Linewidth',8,'Color','k')
line([0 0.5],[0.5 0],'Linewidth',8,'Color','k')
axis off
saveas(gcf,'redT','epsc2');
saveas(gcf,'redT','fig');
