function p = drawGrid(p,lvl)
% draw the mesh

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
myColor = loadField('p.params.output','myColor',p,'k');
drawInfo = loadField('p.params.output','drawInfo',p,true);
holdIt = loadField('p.params.output','holdIt',p,true);

% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

% load enumerated data
nrNodes = p.level(lvl).nrNodes;
nrElems = p.level(lvl).nrElems;

%% draw routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrNodes4e = size(n4e,2);
coordX = reshape(c4n(n4e',1),nrNodes4e,nrElems);
coordY = reshape(c4n(n4e',2),nrNodes4e,nrElems);
val = zeros(size(coordX));

if holdIt
    hold all; 
end

patch(coordX,coordY,val,'FaceColor','none','LineWidth',lineWidth,'EdgeColor',myColor);

hold off

%% draw label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(drawInfo)
	xlabel(sprintf('Nr of Nodes %g',nrNodes));
	title(sprintf('Grid on Level %g',lvl));
end
