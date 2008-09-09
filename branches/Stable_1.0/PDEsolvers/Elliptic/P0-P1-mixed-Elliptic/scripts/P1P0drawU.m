function p = drawU(p,lvl)

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

clf;
hold on
coordX = c4n(:,1);
coordY = c4n(:,2);
Z = [u;u;u];
X = coordX(n4e)';
Y = coordY(n4e)';
fill3(X,Y,Z,Z);

if(drawWalls)
	for curElem = 1:nrElems
		curNodes = n4e(curElem,:);
		coordX = c4n(curNodes([1 2]),1);
		coordY = c4n(curNodes([1 2]),2);

		curU = u(curElem);
		Z = [0,0,curU,curU];
		X = [coordX;coordX([2 1])];
		Y = [coordY;coordY([2 1])];
		fill3(X,Y,Z,[0,0,curU,curU]');

		coordX = c4n(curNodes([2 3]),1);
		coordY = c4n(curNodes([2 3]),2);

		curU = u(curElem);
		Z = [0,0,curU,curU];
		X = [coordX;coordX([2 1])];
		Y = [coordY;coordY([2 1])];
		fill3(X,Y,Z,[0,0,curU,curU]');

		coordX = c4n(curNodes([3 1]),1);
		coordY = c4n(curNodes([3 1]),2);

		curU = u(curElem);
		Z = [0,0,curU,curU];
		X = [coordX;coordX([2 1])];
		Y = [coordY;coordY([2 1])];
		fill3(X,Y,Z,[0,0,curU,curU]');
	end
end
hold off
view(20,40)

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Discrete P0 Solution');
	grid on;
end
