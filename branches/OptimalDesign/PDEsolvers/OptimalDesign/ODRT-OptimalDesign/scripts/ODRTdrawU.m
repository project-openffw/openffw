function p = drawU(p,lvl,drawInfo,drawWalls)
% Copyright 2007 David Guenther
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
%
%
%

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

if(nargin < 3 || isempty(drawInfo))
	drawInfo = true;
end

if(nargin < 4 || isempty(drawWalls))
	drawWalls = true;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

nrElems = p.level(lvl).nrElems;
nrDoF = p.level(lvl).nrDoF;
u = p.level(lvl).u;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cla
hold on
coordX = c4n(:,1);
coordY = c4n(:,2);
Z = [u;u;u];
X = coordX(n4e)';
Y = coordY(n4e)';
fill3(X,Y,Z,Z);

if(drawWalls)
	for curElem = 1:nrElems
		curNodes = n4e(curElem,:);
		coordX = c4n(curNodes([1 2]),1);
		coordY = c4n(curNodes([1 2]),2);

		curU = u(curElem);
		Z = [0,0,curU,curU];
		X = [coordX;coordX([2 1])];
		Y = [coordY;coordY([2 1])];
		fill3(X,Y,Z,[0,0,curU,curU]');

		coordX = c4n(curNodes([2 3]),1);
		coordY = c4n(curNodes([2 3]),2);

		curU = u(curElem);
		Z = [0,0,curU,curU];
		X = [coordX;coordX([2 1])];
		Y = [coordY;coordY([2 1])];
		fill3(X,Y,Z,[0,0,curU,curU]');

		coordX = c4n(curNodes([3 1]),1);
		coordY = c4n(curNodes([3 1]),2);

		curU = u(curElem);
		Z = [0,0,curU,curU];
		X = [coordX;coordX([2 1])];
		Y = [coordY;coordY([2 1])];
		fill3(X,Y,Z,[0,0,curU,curU]');
	end
end
hold off
view(30,30)

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Discrete P0 Solution');
	grid on;
end
