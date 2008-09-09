function p = AWdrawSigma(p,lvl)

if(nargin < 2 || isempty(lvl))
	lvl = p.level(end).level;
end

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set graphic options from structure p
drawInfo = loadField('p.params.output','drawInfo',p,true);
myColor = loadField('p.params.output','myColor',p,'k');
lineWidth = loadField('p.params.output','lineWidth',p,1);
factor = loadField('p.params.output','factor',p,1000);
showGrid = loadField('p.params.output','showGrid',p,true);

% load geometry
c4n = p.level(lvl).geom.c4n;
n4e = p.level(lvl).geom.n4e;

% load discrete solution
dev = p.level(end).dev4n;
u = p.level(end).u4e;

% load enumerated data
nrElems = p.level(lvl).nrElems;
nrDoF = p.level(lvl).nrDoF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uU = squeeze(u(:,1,:));
uV = squeeze(u(:,2,:));
coordX = reshape(c4n(n4e',1),3,nrElems) + factor*uU;
coordY = reshape(c4n(n4e',2),3,nrElems) + factor*uV;
Z = -reshape(dev(n4e'),3,nrElems);

if showGrid
    patch(coordX,coordY,Z,'linewidth',lineWidth,'EdgeColor',myColor);
else
    patch(coordX,coordY,Z);
end

axis equal

if(drawInfo)
	xlabel(sprintf('Nr of degrees of freedom: %g',nrDoF));
	title('Displaced grid and the deviatoric part of the strain');
end
