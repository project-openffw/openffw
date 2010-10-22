function p = ODP2drawU(p,lvl,drawInfo)

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

if(nargin < 3 || isempty(drawInfo))
	drawInfo = true;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = p.level(lvl).u;
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
nrDoF = p.level(lvl).nrDoF;
nrNodes = p.level(lvl).nrNodes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = u(1:nrNodes);
trisurf(n4e,c4n(:,1),c4n(:,2),u','facecolor','interp');

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Discrete P2 Solution');
	grid on;
end
