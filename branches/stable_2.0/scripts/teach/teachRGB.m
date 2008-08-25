function teachRGB
% shows red, green and blue refinement rules

% Copyright 2007 Jan Reininghaus, David Guenther
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


lineWidth = 2;

c4n = [-1 0; 1 0; 0 1; 0 0; 1/2 1/2; -1/2 1/2];


% plot triangle
figure
hold on
n4e = [1 2 3];
[xPatch,yPatch,xRef,yRef]=coord4Plot(c4n,n4e);
patch(xPatch',yPatch','w');
plot(xRef',yRef','k','LineWidth',lineWidth);
title('initial triangle')
hold off

% plot green
figure
hold on

n4e = [3 1 4 ; 2 3 4];
[xPatch,yPatch,xRef,yRef]=coord4Plot(c4n,n4e);
patch(xPatch',yPatch','g');
plot(xRef',yRef','k','LineWidth',lineWidth);
title('green refinement')

hold off


% plot blueR
figure
hold on

n4e=[4 2 5;3 1 4;3 4 5];
[xPatch,yPatch,xRef,yRef]=coord4Plot(c4n,n4e);
patch(xPatch',yPatch','b');
plot(xRef',yRef','k','LineWidth',lineWidth);
title('blue refinement')

hold off


% plot red
figure
hold on

n4e=[1 4 6;5 6 4;4 2 5;6 5 3];
[xPatch,yPatch,xRef,yRef]=coord4Plot(c4n,n4e);
patch(xPatch',yPatch','r');
plot(xRef',yRef','k','LineWidth',lineWidth);
title('red refinement')

hold off

function [xPatch,yPatch,refx,refy]=coord4Plot(c4n,n4e)
nrElems=size(n4e,1);

cx=c4n(:,1);
cy=c4n(:,2);
xPatch = cx(n4e);
yPatch = cy(n4e);
refx=[];
refy=[];
for j = 1:nrElems
	p=c4n(n4e(j,[1,2]),:);
	p1=p(1,:);
	p2=p(2,:);
	t=p2-p1;
	v=[0 -1; 1 0]*t';
	refx(j,:)=c4n(n4e(j,[1,2]),1)+[1 1;1 -1]*[0.08*v(1);0.2*t(1)];
	refy(j,:)=c4n(n4e(j,[1,2]),2)+[1 1;1 -1]*[0.08*v(2);0.2*t(2)];
end

