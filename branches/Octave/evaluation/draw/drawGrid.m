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

drawNodeNr = loadField('p.params.output','drawNodeNr',p,false);
drawElemNr = loadField('p.params.output','drawElemNr',p,false);
drawEdgeNr = loadField('p.params.output','drawEdgeNr',p,false);

% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

% load enumerated data
nrNodes = p.level(lvl).nrNodes;
nrElems = p.level(lvl).nrElems;
nrEdges = p.level(lvl).nrEdges;
midPoint4e = p.level(lvl).enum.midPoint4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;

%% draw routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(n4e,2) == 9
    nrNodes4e = 4;
else    
    nrNodes4e = size(n4e,2);
end

coordX = reshape(c4n(n4e(:,1:nrNodes4e)',1),nrNodes4e,nrElems);
coordY = reshape(c4n(n4e(:,1:nrNodes4e)',2),nrNodes4e,nrElems);
val = zeros(size(coordX));

if holdIt
    hold all; 
end

patch(coordX,coordY,val,'FaceColor','none','LineWidth',lineWidth,'EdgeColor',myColor);


if drawNodeNr
    text(c4n(:,1),c4n(:,2),num2str((1:nrNodes)'), ...
        'HorizontalAlignment','center', ...
        'EdgeColor','k', ...
        'BackgroundColor','r', ...
        'LineStyle','-', ...
        'LineWidth',1, ...
        'Margin',1);
end

if drawElemNr
    text(midPoint4e(:,1),midPoint4e(:,2),num2str((1:nrElems)'), ...
        'HorizontalAlignment','center', ...
        'EdgeColor','k', ...
        'BackgroundColor','b', ...
        'LineStyle','-', ...
        'LineWidth',1, ...
        'Margin',1);
end

if drawEdgeNr
    text(midPoint4ed(:,1),midPoint4ed(:,2),num2str((1:nrEdges)'), ...
        'HorizontalAlignment','center', ...
        'EdgeColor','k', ...
        'BackgroundColor','y', ...
        'LineStyle','-', ...
        'LineWidth',1, ...
        'Margin',1);
end

hold off



% nrNodes4e = size(n4e,2);
% coordX = reshape(c4n(n4e',1),nrNodes4e,nrElems);
% coordY = reshape(c4n(n4e',2),nrNodes4e,nrElems);
% val = zeros(size(coordX));

%if holdIt
%    hold on;
%    hold all; 
%end

%triplot2(n4e,c4n(:,1),c4n(:,2));

%hold off

%% draw label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(drawInfo)
	xlabel(sprintf('Nr of Nodes %g',nrNodes));
	title(sprintf('Grid on Level %g',lvl));
end
