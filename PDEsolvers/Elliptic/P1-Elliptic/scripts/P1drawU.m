function p = P1drawU(p,lvl)
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

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set graphic options from structure p
drawInfo = loadField('p.params.output','drawInfo',p,true);

% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
% n4ed = p.level(lvl).enum.n4ed;

% load discrete solution
u = p.level(lvl).u;

% load enumerated data
nrDoF = p.level(lvl).nrDoF;
nrElems = p.level(lvl).nrElems;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = reshape(u(n4e'),[3 nrElems]);
coordX = reshape(c4n(n4e',1),3,nrElems);
coordY = reshape(c4n(n4e',2),3,nrElems);
patch(coordX,coordY,u);
% triplot3(n4ed,c4n(:,1),c4n(:,2),u);

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title(strcat('Discrete P1 Solution on Level ',num2str(lvl)));
	grid on;
end
