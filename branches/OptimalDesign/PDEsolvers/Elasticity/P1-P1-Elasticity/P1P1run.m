function p = P1P1run(p)
%run.m makes available all necessary initial data,
%handels the discrete displacement,stress, ..., and starts
%the computation for a P1-FE method in linear elasticity. 
%
%authors: David Guenther, Jan Reininghaus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial enumeration
enumerate = p.statics.enumerate;
p = enumerate(p);

p.PDE.mu = p.PDE.E/(2*(1+p.PDE.nu));
p.PDE.lambda = p.PDE.E * p.PDE.nu /( (1+p.PDE.nu)*(1-2*p.PDE.nu) );

% Set initial values
nrNodes = p.level(1).nrNodes;
p.level(1).x = zeros(2*nrNodes,2);

postProc = p.statics.postProc;
p = postProc(p);

p.statics.basisU = @getDisplacementBasis;
p.statics.u_h = @getU_h;
p.statics.sigma_h = @getSigma_h;

p = computeSolution(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_h = getU_h(x,y,curElem,lvl,p)

u = p.level(lvl).u4e;
basisU = p.statics.basisU;
evalBasisU = basisU(x,y,curElem,lvl,p);
nrBasisFuncU = size(evalBasisU,3);
curU = u(:,:,curElem);
curU = curU(:)';
curU = repmat(curU,2*length(x),1);
curU = reshape(curU,[length(x),2,nrBasisFuncU]);

u_h = evalBasisU.*curU;
u_hU = u_h(:,:,1:nrBasisFuncU/2);
u_hV = u_h(:,:,nrBasisFuncU/2+1:end);

u_hU = squeeze(sum(u_hU,3));
u_hV = squeeze(sum(u_hV,3));

u_h = [u_hU(:,1)'; u_hV(:,2)'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigma_h = getSigma_h(x,y,curElem,lvl,p)

sigma = p.level(lvl).sigma;
curSigma = sigma(:,:,curElem);
sigma_h = repmat(curSigma,[1,1,length(x)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basisU = getDisplacementBasis(x,y,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
area4e = p.level(lvl).enum.area4e;

nodes = n4e(curElem,:);
coords = c4n(nodes,:);
area = area4e(curElem);

P1 = coords(1,:);
P2 = coords(2,:);
P3 = coords(3,:);

x = x';
y = y';

b1 = 1/2/area*( (P2(2)-P3(2))*x + (P3(1)-P2(1))*y + P2(1)*P3(2)-P3(1)*P2(2) );
b2 = 1/2/area*( (P3(2)-P1(2))*x + (P1(1)-P3(1))*y + P3(1)*P1(2)-P1(1)*P3(2) );
b3 = 1/2/area*( (P1(2)-P2(2))*x + (P2(1)-P1(1))*y + P1(1)*P2(2)-P2(1)*P1(2) );

basisU = zeros(length(x),2,6);
basisU(:,1,1) = b1;
basisU(:,1,2) = b2;
basisU(:,1,3) = b3;
basisU(:,2,4) = b1;
basisU(:,2,5) = b2;
basisU(:,2,6) = b3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
