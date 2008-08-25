function p = drawMarkedGrid(p,lvl)
% draw all marked edges

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

% set graphic options from structure p
lineWidth = loadField('p.params.output','lineWidth',p,0.1);
myColor = loadField('p.params.output','myColor',p,'m');
drawInfo = loadField('p.params.output','drawInfo',p,true);
holdIt = loadField('p.params.output','holdIt',p,true);

% load geometry
c4n = p.level(lvl).geom.c4n;

% load enumerated data
nrNodes = p.level(lvl).nrNodes;
n4ed = p.level(lvl).enum.n4ed;
markedEdges = p.level(lvl).markedEdges;


%% draw routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~holdIt 
    clf
    drawGrid(p,lvl); 
end

hold on

for j = find(markedEdges)
	plot(c4n(n4ed(j,:),1),c4n(n4ed(j,:),2),myColor,'LineWidth',lineWidth);
end

hold off

%% draw label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(drawInfo)
	xlabel(sprintf('Nr of Nodes %g',nrNodes));
	title(sprintf('Marked Edges after Closure (Level = %g)',lvl));
end
