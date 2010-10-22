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


normRes = localRes*10;

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

nrElems = p.level(lvl).nrElems;
area4e = p.level(lvl).enum.area4e;

nrDoF = p.level(lvl).nrDoF;
% grad = p.level(lvl).grad;
uh = p.statics.u_h;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dofU4e = p.level(lvl).enum.dofU4e;
% x = [1 0 1 1 -1 0 0 0 0 0 0];
% x = rand(1,max(max(dofU4e)));
% x = zeros(1,max(max(dofU4e)));
% x(12) = 1;
% x(20) = 1;
% x(30) = 1;
% x(25) = 1;
% x(15) = 1;

% p.level(lvl).u4e =  x(dofU4e);

u4e = p.level(lvl).u4e;

if(nrElems > 20)
	warning('This may take forever. Please try a smaller level')
	pause
end

clf;
hold on;
for curElem = 1:nrElems
% curElem = 1;
% 	curU = grad(:,1,curElem);
% 	curV = grad(:,2,curElem);
   
	curArea = area4e(curElem);
	curNodes = n4e(curElem,:);
	coordX = c4n(curNodes,1);
	coordY = c4n(curNodes,2);
	minX = min(coordX);
	minY = min(coordY);
	maxX = max(coordX);
	maxY = max(coordY);
    
    node1 = uh(coordX(1),coordY(1),curElem,lvl,p);
    node2 = uh(coordX(2),coordY(2),curElem,lvl,p);
    node3 = uh(coordX(3),coordY(3),curElem,lvl,p);
    curU = [node1(1);node2(1);node3(1)];
	curV = [node1(2);node2(2);node3(2)];
% 	C = sqrt(curU.^2+curV.^2);
% 	fill(coordX,coordY,C);	

	
	diam = sqrt(curArea);
	XI = minX:diam/localRes:maxX;
	YI = minY:diam/localRes:maxY;
	
	UI=tri2grid([coordX,coordY]',[1,2,3,1]' ,...
				curU,XI,YI' ); 
	VI=tri2grid([coordX,coordY]',[1,2,3,1]' ,...
				curV,XI,YI' ); 
			
	XInorm = minX:diam/normRes:maxX;
	YInorm = minY:diam/normRes:maxY;
	
	UInorm=tri2grid([coordX,coordY]',[1,2,3,1]' ,...
				curU,XInorm,YInorm' ); 
	VInorm=tri2grid([coordX,coordY]',[1,2,3,1]' ,...
				curV,XInorm,YInorm' ); 
			
				
				
	C = sqrt(UInorm.^2+VInorm.^2);
	surf(XInorm,YInorm,zeros(size(C)),C,'EdgeColor','none');
					
	% Draw Streamlines
	handle = streamslice(XI,YI,UI,VI,'cubic');
	set(handle,'Color','black');
end
shading interp
% Plot Grid
% drawGrid(p,lvl,'black',2,drawInfo);

hold off

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Gradient of P0-RT0 Solution');
end

