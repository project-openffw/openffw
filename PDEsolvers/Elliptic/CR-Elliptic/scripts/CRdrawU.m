function p = CRdrawU(p,lvl)
% draw displacement

% Copyright 2007 Jan Reininghaus, David Guenther, Andreas Byfut
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

% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

x=c4n(:,1);
y=c4n(:,2);

% load discrete solution
u = p.level(lvl).x;

% load enumerated data
ed4e = p.level(lvl).enum.ed4e;
nrElems = p.level(lvl).nrElems;
nrDoF = p.level(lvl).nrDoF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z = zeros(nrElems,3);
genericT = [-1 1 1; 1 -1 1; 1 1 -1];
for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curU = u(curEdges);
	newU = genericT*curU;
	Z(curElem,:) = newU([2 3 1]);
end

keyboard

coordX = reshape(c4n(n4e',1),3,nrElems);
coordY = reshape(c4n(n4e',2),3,nrElems);
patch(coordX,coordY,Z');

% hold on; 
% for i = 1 : size(n4e,1)
%     curElem =n4e(i,[1 2 3 1]);
%     plot3(x(curElem),y(curElem),zeros(size(x(curElem))),'b');
%     plot3(x(curElem),y(curElem),Z(i,[1 2 3 1]),'r');
% end
% hold off

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Discrete CR Solution');
	grid on;
end


