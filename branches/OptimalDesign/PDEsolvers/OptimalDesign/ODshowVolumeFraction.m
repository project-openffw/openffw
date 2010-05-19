function p = showVolumeFraction(p,lvl,font4axes)
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

if nargin < 3
    font4axes = false;
end

volumeFraction = p.statics.volumeFraction;
nrElems = p.level(lvl).nrElems;
midPoint4e = p.level(lvl).enum.midPoint4e;
c4n = p.level(lvl).geom.c4n;
n4e = p.level(lvl).geom.n4e;
area4e = p.level(lvl).enum.area4e;

clf
hold on

normRes = 13;

coordsX = reshape(c4n(n4e,1),[],3)';
coordsY = reshape(c4n(n4e,2),[],3)';
% 
% for curElem = 1:nrElems;
% 	curArea = area4e(curElem);
% 	curNodes = n4e(curElem,:);
% 	coordX = c4n(curNodes,1);
% 	coordY = c4n(curNodes,2);
% 	minX = min(coordX);
% 	minY = min(coordY);
% 	maxX = max(coordX);
% 	maxY = max(coordY);
% 
%   	diam = sqrt(curArea);
% 
%     XI = minX:diam/normRes:maxX;
% 	YI = minY:diam/normRes:maxY;
% 			
% 	C = volumeFraction(coordsX(:,curElem),coordsY(:,curElem),curElem,lvl,p);
%     
%     C = tri2grid( [coordX,coordY]',[1,2,3,1]',C,XI,YI' );             
% 	surf(XI,YI,zeros(size(C)),C,'EdgeColor','none');
% end

val = zeros(3,nrElems);

for curElem = 1:nrElems
    points = [coordsX(:,curElem), coordsY(:,curElem)];
    val(:,curElem) = volumeFraction(points,curElem,lvl,p);
end

h = patch(coordsX,coordsY,val,val,'FaceColor','interp','EdgeColor','interp');
% set(h,'FaceLighting','phong');


if font4axes
    h = gca;
    set(h,'FontSize',18);
end
