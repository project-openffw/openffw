function p = showVolumeFraction(p,lvl,font4axes)

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
    val(:,curElem) = volumeFraction(coordsX(:,curElem),coordsY(:,curElem),curElem,lvl,p);
end

h = patch(coordsX,coordsY,val,val,'FaceColor','interp','EdgeColor','interp');
% set(h,'FaceLighting','phong');


if font4axes
    h = gca;
    set(h,'FontSize',18);
end