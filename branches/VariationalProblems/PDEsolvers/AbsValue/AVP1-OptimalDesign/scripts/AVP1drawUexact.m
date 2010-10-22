function p = drawUexact(p,lvl,drawInfo)

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

if(nargin < 3 || isempty(drawInfo))
	drawInfo = true;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = p.problem.u_exact;
n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
nrDoF = p.level(lvl).nrDoF;
nrDoF = p.level(lvl).nrDoF;
nrElems = p.level(lvl).nrElems;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = u(c4n(:,1),c4n(:,2));
trisurf(n4e,c4n(:,1),c4n(:,2),u','facecolor','interp');

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Exact P1 Solution');
	grid on;
end


%u = reshape(u(n4e'),[3 nrElems]);
%coordX = reshape(c4n(n4e',1),3,nrElems);
%coordY = reshape(c4n(n4e',2),3,nrElems);
%patch(coordX,coordY,u,u);

%if(drawInfo)
%	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
%	title('Exact P1 Solution');
%	grid on;
%end
