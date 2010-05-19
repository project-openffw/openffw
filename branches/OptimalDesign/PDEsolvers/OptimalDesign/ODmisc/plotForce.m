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


figure(1)
clf
hold on

plot([0 1 1 0 0],[0 0 1 1 0],'k--','LineWidth',2)
plot([0 0.3],[1 1.5],'k--','LineWidth',2)
plot([1 1.3],[1 1.5],'k--','LineWidth',2)
plot([1 1.3],[0 0.5],'k--','LineWidth',2)
plot([0.3 1.3],[1.5 1.5],'k--','LineWidth',2)
plot([1.3 1.3],[0.5 1.5],'k--','LineWidth',2)

plot([0 1 1.3 0.3 0],[0 0 1 1 0],'k-','LineWidth',3)
plot([0.3 0.6],[1 1.5],'k-','LineWidth',3)
plot([1.3 1.6],[1 1.5],'k-','LineWidth',3)
plot([0.6 1.6],[1.5 1.5],'k-','LineWidth',3)
plot([0.6 1.6],[1.5 1.5],'k-','LineWidth',3)
plot([1 1.3],[0 0.5],'k-','LineWidth',3)
plot([1.3 1.6],[0.5 1.5],'k-','LineWidth',3)

text(0.5,1.2,'A','Fontsize',30)
% quiver([1,2],[1.3 1.3],[1.3 1.3],[0 0],'Color','k','Linewidth',3)
annotation(gcf,'arrow',[0.6 0.9],[0.7479 0.7476],'Color','k','Linewidth',4);
text(1.2,1.25,'F','Fontsize',30)

xx = [0:0.01:0.15];
val = spline([0 0.07 0.15],[0.5 0.55 0.5],xx);
plot(xx,val,'Color','k','Linewidth',3);
text(0.06,0.6,'\alpha','Fontsize',30)

axis off
axis auto
