function p = drawRGB(p,lvl)
% draw each triangle in the color how it will be refined
% ( red, green, blue or nothing)

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

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
ed4e = p.level(lvl).enum.ed4e;
newNode4ed = p.level(lvl).enum.newNode4ed;
markedEdges = p.level(lvl).markedEdges;
nrElems = p.level(lvl).nrElems;

%% draw routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newNode4e = newNode4ed(ed4e);
markedElems = find( any(newNode4e,2) );
nrMarkedEd4MarkedElems = sum(markedEdges( ed4e(markedElems,:) ),2);

I = find(nrMarkedEd4MarkedElems == 1);
gElems = markedElems(I);
I = find(nrMarkedEd4MarkedElems == 2);
bElems = markedElems(I);
I = find(nrMarkedEd4MarkedElems == 3);
rElems = markedElems(I);

nrNodes4e = size(n4e,2);
coordX = reshape(c4n(n4e',1),nrNodes4e,nrElems);
coordY = reshape(c4n(n4e',2),nrNodes4e,nrElems);
C = ones([3 nrElems 3]);
C(:,rElems,1) = 1;C(:,rElems,2) = 0;C(:,rElems,3) = 0;
C(:,gElems,1) = 0;C(:,gElems,2) = 1;C(:,gElems,3) = 0;
C(:,bElems,1) = 0;C(:,bElems,2) = 0;C(:,bElems,3) = 1;

if ~holdIt 
    clf
end

hold on

patch(coordX(:,rElems),coordY(:,rElems),'r');
patch(coordX(:,gElems),coordY(:,gElems),'g');
patch(coordX(:,bElems),coordY(:,bElems),'b');

drawRefEdges(p,lvl); 
drawMarkedGrid(p,lvl); 
hold off

if ~holdIt
    p.params.output.holdIt = true;
    drawGrid(p,lvl);
    drawMarkedGrid(p,lvl); 
    drawRefEdges(p,lvl); 
    p.params.output.holdIt = false;
end

view(0,90)

%% title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title(sprintf('RGB (Level = %g)',lvl));

