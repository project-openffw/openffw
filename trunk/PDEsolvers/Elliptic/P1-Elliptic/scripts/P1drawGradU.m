function p = P1drawGradU(p,lvl)
% draw piecewise constant vector field with default parameters:
% drawInfo = true
% localRes = 5
% drawStreamslice = true

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

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set graphic options from structure p
drawInfo = loadField('p.params.output.drawGradU','drawInfo',p,true);
localRes = loadField('p.params.output.drawGradU','localRes',p,5);
drawStreamslice = loadField('p.params.output.drawGradU','drawStreamslice',p,true);

% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

% load gradients
grad = p.level(lvl).gradU4e;

% load enumerated data
nrElems = p.level(lvl).nrElems;
nrDoF = p.level(lvl).nrDoF;
area4e = p.level(lvl).enum.area4e;

%% drawGradU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nrElems > 50)
    warning('This may take forever. Please try a smaller level')
    pause
end

clf;
hold on;

coordX = reshape(c4n(n4e',1),3,nrElems);
coordY = reshape(c4n(n4e',2),3,nrElems);
normGrad = (grad(:,1).^2 + grad(:,2).^2).^(1/2);
patch(coordX,coordY,zeros(size(coordX)),(normGrad*[1 1 1])');

if(drawStreamslice)
    for curElem = 1:nrElems
        curU = grad(curElem,1);
        curV = grad(curElem,2);
        curArea = area4e(curElem);
        curNodes = n4e(curElem,:);
        coordX = c4n(curNodes,1);
        coordY = c4n(curNodes,2);
        minX = min(coordX);
        minY = min(coordY);
        maxX = max(coordX);
        maxY = max(coordY);

        diam = sqrt(curArea);
        XI = minX:diam/localRes:maxX;
        YI = minY:diam/localRes:maxY;

        UI=tri2grid([coordX,coordY]',[1,2,3,1]' ,...
            [curU;curU;curU],XI,YI' );
        VI=tri2grid([coordX,coordY]',[1,2,3,1]' ,...
            [curV;curV;curV],XI,YI' );

        % Draw Streamlines
        handle = streamslice(XI,YI,UI,VI);
        set(handle,'Color','black');
    end
end
hold off

if(drawInfo)
    xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
    title('Piecewise constant gradient of P1 Solution');
end