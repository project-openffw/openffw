function p = AVP1drawUminusUh(p,lvl,drawInfo)

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

if(nargin < 3 || isempty(drawInfo))
	drawInfo = true;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = p.level(lvl).u;
u_exact = p.problem.u_exact;
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
nrDoF = p.level(lvl).nrDoF;
nrNodes = p.level(lvl).nrNodes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_diff = abs(u(1:nrNodes)-u_exact(c4n(:,1),c4n(:,2)));
trisurf(n4e,c4n(:,1),c4n(:,2),u_diff','facecolor','interp');

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Difference Exact and P1 Solution');
	grid on;
end
