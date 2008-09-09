function p = CRrun(p)
%run.m makes available all necessary initial data,
%handels the discrete displacement,flux, ..., and starts
%the computation for a nonconforming CR-FE method. 
%
%authors: David Guenther, Jan Reininghaus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial enumeration
enumerate = p.statics.enumerate;
p = enumerate(p);

% Set initial values
nrEdges = p.level(1).nrEdges;
p.level(1).x = zeros(nrEdges,1);

p.statics.basisU = @getDisplacementBasis;
p.statics.u_h = @getU_h;
p.statics.sigma_h = @getSigma_h;

p = CRpostProc(p);
p = computeSolution(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_h = getU_h(x,y,curElem,lvl,p)

u = p.level(lvl).u4e;
basisU = p.statics.basisU;
basisU = basisU(x,y,curElem,lvl,p);

curU = u(curElem,:);
u_h = curU * basisU';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigma_h = getSigma_h(x,y,curElem,curKappa,lvl,p)

grad4e = p.level(lvl).grad4e;
curGrad = grad4e(curElem,:);

curGrad = reshape(curGrad'*ones(1,length(x)),[2 1 length(x)]);
sigma_h = matMul(curKappa,curGrad);
sigma_h = reshape(sigma_h,2,[])';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basisU = getDisplacementBasis(x,y,curElem,lvl,p)

ed4e = p.level(lvl).enum.ed4e;
midPoint4ed = p.level(lvl).enum.midPoint4ed;
area4e = p.level(lvl).enum.area4e;

edges = ed4e(curElem,:);
coords = midPoint4ed(edges,:);
area = area4e(curElem);

area = area/4;

P1 = coords(1,:);
P2 = coords(2,:);
P3 = coords(3,:);

x = x';
y = y';

b1 = 1/2/area*( (P2(2)-P3(2))*x + (P3(1)-P2(1))*y + P2(1)*P3(2)-P3(1)*P2(2) );
b2 = 1/2/area*( (P3(2)-P1(2))*x + (P1(1)-P3(1))*y + P3(1)*P1(2)-P1(1)*P3(2) );
b3 = 1/2/area*( (P1(2)-P2(2))*x + (P2(1)-P1(1))*y + P1(1)*P2(2)-P2(1)*P1(2) );

basisU = [b1;b2;b3]';

