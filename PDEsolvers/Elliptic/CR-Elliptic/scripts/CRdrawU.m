function p = CRdrawU(p,lvl)

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set graphic options from structure p
drawInfo = loadField('p.params.output','drawInfo',p,true);

% load geometry
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;

% load discrete solution
u = p.level(lvl).x;

% load enumerated data
ed4e = p.level(lvl).enum.ed4e;
nrElems = p.level(lvl).nrElems;
nrDoF = p.level(lvl).nrDoF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z = zeros(3,nrElems);
genericT = [-1 1 1; 1 -1 1; 1 1 -1];
for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curU = u(curEdges);
	newU = genericT*curU;
	Z(:,curElem) = newU([2 3 1]);
end

coordX = reshape(c4n(n4e',1),3,nrElems);
coordY = reshape(c4n(n4e',2),3,nrElems);
patch(coordX,coordY,Z,Z);

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Discrete CR Solution');
	grid on;
end


