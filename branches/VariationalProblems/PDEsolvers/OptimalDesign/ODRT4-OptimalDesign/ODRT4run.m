function p = ODRTrun(p)
% author: David Guenther 

%% initialization
% Initial enumeration
p = ODRT4enumerate(p);

% Set initial values
nrEdges = p.level(1).nrEdges;
nrElems = p.level(1).nrElems;
p.level(1).x = zeros(nrEdges + nrElems,1);

% supply the basis functions
p.statics.basisU = @getDisplacementBasis;
p.statics.stressBasis = @getStressBasis;
p.statics.divStressBasis = @getDivStressBasis;
p.statics.jacobianStressBasis = @getJacobianStressBasis;
p.statics.basisP1 = @getP1basis;
p.statics.gradBasisP1 = @getGradP1basis;
p.statics.basisP1NC = @getP1NCbasis;
p.statics.gradBasisP1NC = @getGradP1NCbasis;
p.statics.basisP2 = @getP2basis;
p.statics.gradBasisP2 = @getGradP2basis;

p = ODRT4getUhFuncs(p);
p = ODRT4getPhFuncs(p);

p.statics.basisCoefficients = ...
    [ 1 0 0 -2  0 -2 
      0 1 0 -2 -2  0
      0 0 1  0 -2 -2
      0 0 0  4  0  0
      0 0 0  0  4  0
      0 0 0  0  0  4 ] ;
  
  
  
% compute the discrete solution
p = ODcomputeSolution(p);


%% supply basis functions for p_{\epsilon,h}
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

signum = ones(3,1);
I = find(e4ed(edges,2) == curElem);
signum(I) = -1;

b1 = signum(1)*length4ed(edges(1))/(2*area)*([x-coords(3,1),y-coords(3,2)] );	
b2 = signum(2)*length4ed(edges(2))/(2*area)*([x-coords(1,1),y-coords(1,2)] );	
b3 = signum(3)*length4ed(edges(3))/(2*area)*([x-coords(2,1),y-coords(2,2)] );

basisSigma = zeros(3,2,length(x));

basisSigma(1,:,:) = b1';
basisSigma(2,:,:) = b2';
basisSigma(3,:,:) = b3';

%% supply jacobian of the basis functions for p_{\epsilon,h}
function jacobianBasisSigma = getJacobianStressBasis(x,y,curElem,lvl,p)

c4n = p.level(lvl).geom.c4n;
n4e = p.level(lvl).geom.n4e;
ed4e = p.level(lvl).enum.ed4e;
e4ed = p.level(lvl).enum.e4ed;
area4e = p.level(lvl).enum.area4e;
length4ed = p.level(lvl).enum.length4ed;

edges = ed4e(curElem,:);
coords = c4n(n4e(curElem,:),:);
area = area4e(curElem);

signum = ones(3,1);
I = find(e4ed(edges,2) == curElem);
signum(I) = -1;

b1x = signum(1)*length4ed(edges(1))/(2*area);	
b2x = signum(2)*length4ed(edges(2))/(2*area);	
b3x = signum(3)*length4ed(edges(3))/(2*area);

b1y = signum(1)*length4ed(edges(1))/(2*area);	
b2y = signum(2)*length4ed(edges(2))/(2*area);	
b3y = signum(3)*length4ed(edges(3))/(2*area);

jacobianBasisSigma = zeros(2,2,3,length(x));

jacobianBasisSigma(:,:,1,:) = repmat([b1x 0;0 b1y],[1 1 length(x)]);
jacobianBasisSigma(:,:,2,:) = repmat([b2x 0;0 b2y],[1 1 length(x)]);
jacobianBasisSigma(:,:,3,:) = repmat([b3x 0;0 b3y],[1 1 length(x)]);

%% supply the divergence of the basis functions of p_{\epsilon,h}
function divBasisSigma = getDivStressBasis(x,y,curElem,lvl,p)

ed4e = p.level(lvl).enum.ed4e;
e4ed = p.level(lvl).enum.e4ed;
area4e = p.level(lvl).enum.area4e;
length4ed = p.level(lvl).enum.length4ed;

edges = ed4e(curElem,:);
area = area4e(curElem);

signum = ones(3,1);
I = find(e4ed(edges,2) == curElem);
signum(I) = -1;

divBasisSigma = zeros(3,1,length(x));

for j = 1:length(x)
    b1 = signum(1)*length4ed(edges(1))/area;	
    b2 = signum(2)*length4ed(edges(2))/area;	
    b3 = signum(3)*length4ed(edges(3))/area;

    divBasisSigma(:,:,j) = [b1';b2';b3'];
end

%% supply the basis function for u_h
function basisU = getDisplacementBasis(x,y,curElem,lvl,p)

basisU = ones(length(x),1);

%% supply the P1-basis functions
function basisP1 = getP1basis(x,y,curElem,lvl,p)

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

basisP1 = [b1;b2;b3]';

%% supply the gradient of the P1-basis functions
function gradBasisP1 = getGradP1basis(x,y,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
area4e = p.level(lvl).enum.area4e;

nodes = n4e(curElem,:);
coords = c4n(nodes,:);
area = area4e(curElem);

P1 = coords(1,:);
P2 = coords(2,:);
P3 = coords(3,:);

b1 = 1/2/area*[ (P2(2)-P3(2)), (P3(1)-P2(1)) ];
b2 = 1/2/area*[ (P3(2)-P1(2)), (P1(1)-P3(1)) ];
b3 = 1/2/area*[ (P1(2)-P2(2)), (P2(1)-P1(1)) ];

gradBasisP1 = [b1;b2;b3];

gradBasisP1 = repmat(gradBasisP1,[1 1 length(x)]);

%% supply P1NC - basis
function basisP1NC = getP1NCbasis(x,y,curElem,lvl,p)

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

basisP1NC = [b1;b2;b3]';

%% supply the gradient of the P1NC-basis functions
function gradBasisP1NC = getGradP1NCbasis(x,y,curElem,lvl,p)

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

b1 = 1/2/area*[ (P2(2)-P3(2)), (P3(1)-P2(1)) ];
b2 = 1/2/area*[ (P3(2)-P1(2)), (P1(1)-P3(1)) ];
b3 = 1/2/area*[ (P1(2)-P2(2)), (P2(1)-P1(1)) ];

gradBasisP1NC = repmat([b1;b2;b3],[1,1,length(x)]);

%% supply P2-basis
function basisP2 = getP2basis(x,y,curElem,lvl,p)

basisP1 = p.statics.basisP1;
evalBasisP1 = basisP1(x,y,curElem,lvl,p);

b1 = evalBasisP1(:,1)';
b2 = evalBasisP1(:,2)';
b3 = evalBasisP1(:,3)';

basisP2 = [b1; b2; b3; b1.*b2; b2.*b3; b3.*b1]';

%% supply grad P2-basis
function gradBasisP2 = getGradP2basis(x,y,curElem,lvl,p)

% gradBasisP1 = p.statics.gradBasisP1;
% basisP1 = p.statics.basisP1;
% 
% evalGradBasisP1 = gradBasisP1(x,y,curElem,lvl,p);
% evalBasisP1 = basisP1(x,y,curElem,lvl,p)';
% 
% evalGradBasisP1X = squeeze(evalGradBasisP1(:,1,:));
% evalGradBasisP1Y = squeeze(evalGradBasisP1(:,2,:));
% 
% var1 = evalGradBasisP1X.*evalBasisP1([2 3 1],:) + evalGradBasisP1X([2 3 1],:).*evalBasisP1;
% var2 = evalGradBasisP1Y.*evalBasisP1([2 3 1],:) + evalGradBasisP1Y([2 3 1],:).*evalBasisP1;
% 
% gradBasisP2 = zeros(6,2,length(x));
% gradBasisP2(1:3,:,:) = evalGradBasisP1;
% gradBasisP2(4:6,1,:) = var1;
% gradBasisP2(4:6,2,:) = var2;

P2Basis = p.statics.basisP2;

P1grad4e = p.level(lvl).enum.grad4e;
C = p.statics.basisCoefficients;

curP1Grad = P1grad4e(:,:,curElem);
curBasisU = P2Basis(x,y,curElem,lvl,p);

basisSigma  = zeros(6,2,length(x));

for i = 1 : length(x)
    curGradP2 = [ curP1Grad ;
        curP1Grad(1,:)*curBasisU(i,2) + curBasisU(i,1)*curP1Grad(2,:) ;
        curP1Grad(2,:)*curBasisU(i,3) + curBasisU(i,2)*curP1Grad(3,:) ;
        curP1Grad(1,:)*curBasisU(i,3) + curBasisU(i,1)*curP1Grad(3,:) ] ;
    basisSigma(:,:,i) = C * curGradP2;
end

gradBasisP2 = basisSigma;