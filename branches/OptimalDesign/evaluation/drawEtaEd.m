function p = drawEtaEd(p,lvl)
% draw the local edge error indicators

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
lineWidth = loadField('p.params.output','lineWidth',p,1.5);
nrColors = loadField('p.params.output','nrColors',p,256);
showColorbar = loadField('p.params.output','showColorbar',p,false);

c4n = p.level(lvl).geom.c4n;
n4ed = p.level(lvl).enum.n4ed;
etaEd = p.level(lvl).etaEd;
nrEdges = p.level(lvl).nrEdges;

%% draw routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map = colormap(jet(nrColors));

etaEd = etaEd./max(etaEd);
etaEd = etaEd *(nrColors);

if ~holdIt 
    clf
    drawGrid(p,lvl); 
end

hold on

for j = 1:nrEdges
	colorIdx = floor(etaEd(j));
	if(colorIdx > 0)
		plot(c4n(n4ed(j,:),1),c4n(n4ed(j,:),2),'LineWidth',lineWidth,'Color',map(colorIdx,:));
	end
end

hold off

%% colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if showColorbar
    colorbar('location','SouthOutside','XTickLabel',{'Low Error','','','High Error',''})
end

%% title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title(sprintf('Estimated Error on Edges (Level = %g)',lvl));
