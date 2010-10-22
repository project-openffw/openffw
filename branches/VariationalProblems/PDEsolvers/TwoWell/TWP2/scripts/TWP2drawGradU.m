function p = drawGradU(p,lvl,localRes,drawInfo)

if(nargin < 2)
	lvl = p.level(end).level;
end

if(nargin < 3)
	localRes = 10;
end

if(nargin < 4)
	drawInfo = true;
end


if(isempty(lvl))
	lvl = p.level(end).level;
end

if(isempty(drawInfo))
	drawInfo = true;
end

if(isempty(localRes))
	localRes = 10;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

nrElems = p.level(lvl).nrElems;
area4e = p.level(lvl).enum.area4e;

nrDoF = p.level(lvl).nrDoF;
grad = p.level(lvl).grad;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nrElems > 50)
	warning('This may take forever. Please try a smaller level')
	pause
end

clf;
hold on;
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

	C = sqrt(curU.^2+curV.^2);
	fill(coordX,coordY,C);			
	
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

% Plot Grid
drawGrid(p,lvl,'black',2);

hold off

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Piecewise constant gradient of P2 Solution');
end

