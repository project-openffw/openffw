function p = drawDivU(p,lvl,drawInfo)

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

if(nargin < 3 || isempty(drawInfo))
	drawInfo = false;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

% load discrete solution

divU_h  = p.statics.divU_h ;
% load enumerated data
ed4e = p.level(lvl).enum.ed4e;
nrElems = p.level(lvl).nrElems;
nrDoF = p.level(lvl).nrDoF;
etaT = p.level(lvl).etaT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c4n = c4n +rand(size(c4n))*1e-3;
Z = zeros(3,nrElems);
for curElem = 1:nrElems
	coords = c4n(n4e(curElem,:),:);
divUh = divU_h(coords(:,1),coords(:,2),curElem,lvl,p);
    val = divUh;
%     val = etaT;
    Z(:,curElem) = sqrt(sum(val.*val));
end







coordX = c4n(:,1);
coordY = c4n(:,2);
X = coordX(n4e)';
Y = coordY(n4e)';
fill3(X,Y,Z,Z);	

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Discrete CR Solution');
	grid on;
end



