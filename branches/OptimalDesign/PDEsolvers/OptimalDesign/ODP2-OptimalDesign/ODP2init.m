function p = ODP2init(p)
% author: David Guenther 

%% initialization
% Initial enumeration
enumerate = p.statics.enumerate;
p = enumerate(p);

% Set initial values
nrNodes = p.level(1).nrNodes;
p.level(1).x = zeros(nrNodes,1);

% supply the discrete functions
p.statics.basisU = @getDisplacementBasis;
p.statics.stressBasis = @getStressBasis;

p = ODP2getUhFuncs(p);
p = ODP2getPhFuncs(p);

p.statics.basisCoefficients = ...
    [ 1 0 0 -2  0 -2 
      0 1 0 -2 -2  0
      0 0 1  0 -2 -2
      0 0 0  4  0  0
      0 0 0  0  4  0
      0 0 0  0  0  4 ] ;
  
% compute the discrete solution
% p = ODcomputeSolution(p);

%% supply displacement basis (P1)
function basisU = getDisplacementBasis(points,curElem,lvl,p)

basisU = P2Basis(points,curElem,lvl,p)*p.statics.basisCoefficients';

function basisU = P2Basis(points,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
area4e = p.level(lvl).enum.area4e;

nodes = n4e(curElem,:);
coords = c4n(nodes,:);
area = area4e(curElem);

P1 = coords(1,:);
P2 = coords(2,:);
P3 = coords(3,:);

x = points(:,1)';
y = points(:,2)';

b1 = 1/2/area*( (P2(2)-P3(2))*x + (P3(1)-P2(1))*y + P2(1)*P3(2)-P3(1)*P2(2) );
b2 = 1/2/area*( (P3(2)-P1(2))*x + (P1(1)-P3(1))*y + P3(1)*P1(2)-P1(1)*P3(2) );
b3 = 1/2/area*( (P1(2)-P2(2))*x + (P2(1)-P1(1))*y + P1(1)*P2(2)-P2(1)*P1(2) );

basisU = [b1;
          b2;
          b3;
          b1.*b2;
          b2.*b3;
          b1.*b3]';


%% supply stress basis
function basisSigma = getStressBasis(points,curElem,lvl,p)

P1grad4e = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;

curP1Grad = P1grad4e(:,:,curElem);
curBasisU = P2Basis(points,curElem,lvl,p);

basisSigma  = zeros(6,2,length(points(:,1)));

for i = 1 : length(points(:,1))
    curGradP2 = [ curP1Grad ;
        curP1Grad(1,:)*curBasisU(i,2) + curBasisU(i,1)*curP1Grad(2,:) ;
        curP1Grad(2,:)*curBasisU(i,3) + curBasisU(i,2)*curP1Grad(3,:) ;
        curP1Grad(1,:)*curBasisU(i,3) + curBasisU(i,1)*curP1Grad(3,:) ] ;
    basisSigma(:,:,i) = C * curGradP2;
end
