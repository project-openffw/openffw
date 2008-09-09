function p = RT0P0drawU(p,lvl)

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set graphic options from structure p
drawInfo = loadField('p.params.output','drawInfo',p,true);
drawWalls = loadField('p.params.output','drawWalls',p,true);

% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

% load discrete solution
u = p.level(lvl).u;

% load enumerated data
nrElems = p.level(lvl).nrElems;
nrDoF = p.level(lvl).nrDoF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cla;
hold on

Z = [u;u;u];
coordX = reshape(c4n(n4e',1),3,nrElems);
coordY = reshape(c4n(n4e',2),3,nrElems);
patch(coordX,coordY,Z,Z);

X = zeros(4,3*nrElems);
Y = zeros(4,3*nrElems);
Z = zeros(4,3*nrElems);

if(drawWalls)
	for curElem = 1:nrElems
		curNodes = n4e(curElem,:);
		coordX = c4n(curNodes([1 2]),1);
		coordY = c4n(curNodes([1 2]),2);

		curU = u(curElem);
		Z(:,3*(curElem-1)+1) = [0,0,curU,curU];
		X(:,3*(curElem-1)+1) = [coordX;coordX([2 1])];
		Y(:,3*(curElem-1)+1) = [coordY;coordY([2 1])];

		coordX = c4n(curNodes([2 3]),1);
		coordY = c4n(curNodes([2 3]),2);

		curU = u(curElem);
		Z(:,3*(curElem-1)+2) = [0,0,curU,curU];
		X(:,3*(curElem-1)+2) = [coordX;coordX([2 1])];
		Y(:,3*(curElem-1)+2) = [coordY;coordY([2 1])];

		coordX = c4n(curNodes([3 1]),1);
		coordY = c4n(curNodes([3 1]),2);

		curU = u(curElem);
		Z(:,3*(curElem-1)+3) = [0,0,curU,curU];
		X(:,3*(curElem-1)+3) = [coordX;coordX([2 1])];
		Y(:,3*(curElem-1)+3) = [coordY;coordY([2 1])];
	end
end

patch(X,Y,Z,Z)
hold off
view(30,30)

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Discrete P0 Solution');
	grid on;
end
