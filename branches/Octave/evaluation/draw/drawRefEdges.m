function p = drawRefEdges(p,lvl)
% draw all reference edges

% Copyright 2007 Jan Reininghaus, David Guenther, 
%                Andreas Byfut, Joscha Gedicke
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


if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

%% input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

holdIt = loadField('p.params.output','holdIt',p,true);
lineWidth = loadField('p.params.output','lineWidth',p,0.1);
myColor = loadField('p.params.output','myColor',p,'k');

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
nrElems = p.level(lvl).nrElems;

%% draw routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~holdIt 
    clf
    drawGrid(p,lvl); 
end

hold on

for j = 1:nrElems
	p=c4n(n4e(j,[1,2]),:);
	p1=p(1,:);
	p2=p(2,:);
	t=p2-p1;
	v=[0 -1; 1 0]*t';
	x=c4n(n4e(j,[1,2]),1)+[1 1;1 -1]*[0.08*v(1);0.2*t(1)];
	y=c4n(n4e(j,[1,2]),2)+[1 1;1 -1]*[0.08*v(2);0.2*t(2)];
	plot(x,y,myColor,'LineWidth',1);
	%plot(c4n(n4e(j,[1,2]),1),c4n(n4e(j,[1,2]),2),myColor,'LineWidth',lineWidth);
end


hold off

%% title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title(sprintf('Reference Edges (Level = %g)',lvl));
