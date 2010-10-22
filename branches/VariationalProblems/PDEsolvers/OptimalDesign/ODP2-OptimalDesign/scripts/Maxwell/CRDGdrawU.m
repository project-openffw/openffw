function p = drawU(p,lvl,drawInfo)

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
% u = p.level(lvl).x;
u_h = p.statics.u_h;
divU_h  = p.statics.divU_h ;
% u_exact = p.problem.u_exact;
% divU_exact = p.problem.divU_exact;
% load enumerated data
ed4e = p.level(lvl).enum.ed4e;
nrElems = p.level(lvl).nrElems;
nrDoF = p.level(lvl).nrDoF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c4n = c4n +rand(size(c4n))*1e-3;
Z = zeros(3,nrElems);
for curElem = 1:nrElems
	coords = c4n(n4e(curElem,:),:);
    normCoords = sum(coords.*coords,2);
    if(min(normCoords)<= 1e-3)
        continue;
    end
	uh = u_h(coords(:,1),coords(:,2),curElem,lvl,p);
%     uExact = u_exact(coords(:,1),coords(:,2),p)';
% divUh = divU_h(coords(:,1),coords(:,2),curElem,lvl,p);
%     val = divUh;
%     val = uh;
    val = uh(:,:);
    
%     Z(:,curElem) = val;
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

view(120,60)

