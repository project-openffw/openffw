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
