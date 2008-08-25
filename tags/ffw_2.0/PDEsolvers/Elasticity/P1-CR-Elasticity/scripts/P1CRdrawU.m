function p = P1CRdrawU(p,lvl)
% draw displacement

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


if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set graphic options from structure p
drawInfo = loadField('p.params.output','drawInfo',p,true);
myColor = loadField('p.params.output','myColor',p,'k');
lineWidth = loadField('p.params.output','lineWidth',p,1);
factor = loadField('p.params.output','factor',p,1000);

% load geometry
c4n = p.level(lvl).geom.c4n;
n4e = p.level(lvl).geom.n4e;

% load discrete solution
u = p.level(lvl).u;

% load enumerated data
nrDoF = p.level(lvl).nrDoF;
nrNodes = p.level(lvl).nrNodes;
nrElems = p.level(lvl).nrElems;
ed4e = p.level(lvl).enum.ed4e;
area4e = p.level(lvl).enum.area4e;
area4n = p.level(lvl).enum.area4n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uU = u(1:nrNodes);
uV = u(nrNodes+1:end);


newUV = zeros(nrNodes,1);
genericT = [-1 1 1; 1 -1 1; 1 1 -1];
for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curU = uV(curEdges);
	curNodes = n4e(curElem,:);
	newU = genericT*curU;
	area = area4e(curElem);
	newUV(curNodes) = newUV(curNodes) + newU([2 3 1])*area;
end


uV = newUV ./ area4n;

triplot(n4e,factor*uU+c4n(:,1), factor*uV+c4n(:,2),...
			'linewidth',lineWidth,'Color',myColor);
			
if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Displaced grid of P1-CR solution');
end




