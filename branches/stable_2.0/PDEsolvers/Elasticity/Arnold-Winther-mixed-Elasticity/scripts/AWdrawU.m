function p = AWdrawU(p,lvl)
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
n4e = p.level(lvl).geom.n4e';

% load discrete solution
u = p.level(lvl).u4e;

% load enumerated data
nrElems = p.level(lvl).nrElems;
angles4e = p.level(lvl).enum.angles4e;
angle4n = p.level(lvl).enum.angle4n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uU = u(:,1,:);
uV = u(:,2,:);

uU = reshape(uU,3,nrElems).* angles4e';
uV = reshape(uV,3,nrElems).* angles4e';

uX4n = accumarray(n4e(:),uU(:))./angle4n;
uY4n = accumarray(n4e(:),uV(:))./angle4n;

nrDoF = p.level(lvl).nrDoF;

triplot(n4e',factor*uX4n+c4n(:,1), factor*uY4n+c4n(:,2),...
			'linewidth',lineWidth,'Color',myColor);
			
if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Displaced grid of Arnold-Winther-mixed solution');
end

