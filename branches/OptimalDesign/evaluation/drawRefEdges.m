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
myColor = loadField('p.params.output','myColor',p,'c');

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
nrElems = p.level(lvl).nrElems;

%% draw routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~holdIt 
    clf
    drawGrid(p,lvl); 
end

hold on

x=c4n(:,1);
y=c4n(:,2);
lambda=.9;
lambda2=.9;
x=(lambda.*x(n4e(:,1:2))+(1-lambda).*x(n4e(:,[3;3])))';
y=(lambda.*y(n4e(:,1:2))+(1-lambda).*y(n4e(:,[3;3])))';
plot(lambda2*x+(1-lambda2)*x([2;1],:),...
    lambda2*y+(1-lambda2)*y([2;1],:),myColor,'LineWidth',lineWidth)

%for j = 1:nrElems
%	plot(c4n(n4e(j,[1,2]),1),c4n(n4e(j,[1,2]),2),myColor,'LineWidth',lineWidth);
%end

hold off

%% title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title(sprintf('Reference Edges (Level = %g)',lvl));
