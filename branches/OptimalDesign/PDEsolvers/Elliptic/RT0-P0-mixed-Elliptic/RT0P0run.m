function p = RT0P0run(p)
%RT0P0run.m makes available all necessary initial data,
%handels the discrete displacement,flux, ..., and starts
%the computation for the mixed RT0-P0-FE method. 
%
%authors: David Guenther, Jan Reininghaus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial enumeration
enumerate = p.statics.enumerate;
p = enumerate(p);

% Set initial values
nrEdges = p.level(1).nrEdges;
nrElems = p.level(1).nrElems;
p.level(1).x = zeros(nrEdges + nrElems,1);

p.statics.basisU = @getDisplacementBasis;
p.statics.u_h = @getU_h;
p.statics.sigma_h = @getSigma_h;

p = RT0P0postProc(p);
p = computeSolution(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_h = getU_h(x,y,curElem,lvl,p)

u = p.level(lvl).u4e;
curU = u(curElem);
u_h = repmat(curU,1,length(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigma_h = getSigma_h(x,y,curElem,curKappa,lvl,p)

basisSigma = getStressBasis(x,y,curElem,lvl,p);

grad4e = p.level(lvl).grad4e;

sigma_h = zeros(length(x),2);

for j = 1:length(x)
    sigma_h(j,:) = curKappa(:,:,j)*(grad4e(curElem,:)*basisSigma(:,:,j))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basisSigma = getStressBasis(x,y,curElem,lvl,p)

c4n = p.level(lvl).geom.c4n;
n4e = p.level(lvl).geom.n4e;
ed4e = p.level(lvl).enum.ed4e;
e4ed = p.level(lvl).enum.e4ed;
area4e = p.level(lvl).enum.area4e;
length4ed = p.level(lvl).enum.length4ed;

edges = ed4e(curElem,:);
coords = c4n(n4e(curElem,:),:);
area = area4e(curElem);

P1 = coords(1,:)';
P2 = coords(2,:)';
P3 = coords(3,:)';

signum = ones(3,1);
I = find(e4ed(edges,2) == curElem);
signum(I) = -1;

basisSigma = zeros(3,2,length(x));

for j = 1:length(x)
    curX = x(j);
    curY = y(j);
    
    b1 = signum(1)*length4ed(edges(1))/(2*area)*([curX;curY] - P3);	
    b2 = signum(2)*length4ed(edges(2))/(2*area)*([curX;curY] - P1);	
    b3 = signum(3)*length4ed(edges(3))/(2*area)*([curX;curY] - P2);

    basisSigma(:,:,j) = [b1';b2';b3'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basisU = getDisplacementBasis(x,y,curElem,lvl,p)

basisU = ones(length(x),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

