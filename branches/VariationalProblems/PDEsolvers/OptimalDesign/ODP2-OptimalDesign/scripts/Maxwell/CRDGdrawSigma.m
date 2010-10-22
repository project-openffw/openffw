function p = drawSigma(p,lvl,factor,lineWidth,myColor,drawInfo)
% drawSigma(p,lvl,factor,lineWidth,myColor,drawInfo)

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

if(nargin < 3 || isempty(factor))
	factor = 1000;
end

if(nargin < 4 || isempty(lineWidth))
	lineWidth = 1;
end

if(nargin < 5 || isempty(myColor))
	myColor = 'k';
end

if(nargin < 6 || isempty(drawInfo))
	drawInfo = true;
end

c4n = p.level(lvl).geom.c4n;
n4e = p.level(lvl).geom.n4e;
nrElems = p.level(lvl).nrElems;
dev = p.level(end).dev4n;
Au = p.level(end).Au;

for curElem = 1:nrElems
    curNodes = n4e(curElem,:);
    curU = Au(curNodes,:);
    
    trisurf([1 2 3],    c4n(n4e(curElem,:),1) + factor*curU(:,1),...
                        c4n(n4e(curElem,:),2) + factor*curU(:,2),...
            -dev(n4e(curElem,:))');
    hold on;
    triplot([1 2 3], c4n(n4e(curElem,:),1) + factor*curU(:,1),...
                     c4n(n4e(curElem,:),2) + factor*curU(:,2),...
                     'linewidth',lineWidth,'Color',myColor);
    hold on;
end

% triplot(n4e,c4n(:,1),c4n(:,2),'linewidth',lineWidth/10,'Color','k')

colormap(gray)
shading interp
view(2)
axis equal
hold off
