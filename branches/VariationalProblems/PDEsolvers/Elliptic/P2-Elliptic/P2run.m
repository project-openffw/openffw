function p = P2run(p)

% author: Joscha Gedicke

% Initial enumeration
enumerate = p.statics.enumerate;
p = enumerate(p);

% Set initial values
nrNodes = p.level(1).nrNodes;
nrEdges = p.level(1).nrEdges;
p.level(1).x = zeros(nrNodes+nrEdges,1);

p.statics.basisU = @getDisplacementBasis;
p.statics.basisGradU = @getGradP2;
p.statics.u_h = @getU_h;
p.statics.gradU_h = @getGradU_h;
p.statics.D2U_h = @getD2U_h;
p.statics.sigma_h = @getSigma_h;

p.statics.basisCoefficients = ...
    [ 1 0 0 -2  0 -2 
      0 1 0 -2 -2  0
      0 0 1  0 -2 -2
      0 0 0  4  0  0
      0 0 0  0  4  0
      0 0 0  0  0  4 ] ;

p = P2postProc(p);
p = computeSolution(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_h = getU_h(x,y,curElem,lvl,p)

u = p.level(lvl).u4e;
basisU = p.statics.basisU;
basisU = basisU(x,y,curElem,lvl,p);

curU = u(curElem,:);
u_h = curU * basisU';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gradU_h = getGradU_h(x,y,curElem,lvl,p)

u = p.level(lvl).u4e;
curU = u(curElem,:);

curGradP2 = getGradP2(x,y,curElem,lvl,p);

for i = 1 : length(x)
    gradU_h(i,:) = curU * curGradP2(:,:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gradP2 = getGradP2(x,y,curElem,lvl,p)

P1grad4e = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;

curP1Grad = P1grad4e(:,:,curElem);
curBasisU = P2Basis(x,y,curElem,lvl,p);

gradP2  = zeros(6,2,length(x));

for i = 1 : length(x)
    curGradP2 = [ curP1Grad ;
        curP1Grad(1,:)*curBasisU(i,2) + curBasisU(i,1)*curP1Grad(2,:) ;
        curP1Grad(2,:)*curBasisU(i,3) + curBasisU(i,2)*curP1Grad(3,:) ;
        curP1Grad(1,:)*curBasisU(i,3) + curBasisU(i,1)*curP1Grad(3,:) ] ;
    gradP2(:,:,i) = C * curGradP2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D2U_h = getD2U_h(x,y,curElem,lvl,p)

u = p.level(lvl).u4e;
P1grad4e = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;

curU = u(curElem,:);
curP1Grad = P1grad4e(:,:,curElem);

D2U_h = zeros(2,2);
D2U_h(1,1) = curU * C * [0;0;0;
       2*curP1Grad(1,1)*curP1Grad(2,1);
       2*curP1Grad(2,1)*curP1Grad(3,1);
       2*curP1Grad(1,1)*curP1Grad(3,1)];
D2U_h(2,2) = curU * C * [0;0;0;
       2*curP1Grad(1,2)*curP1Grad(2,2);
       2*curP1Grad(2,2)*curP1Grad(3,2);
       2*curP1Grad(1,2)*curP1Grad(3,2)];
D2U_h(1,2) = curU * C * [0;0;0;
       curP1Grad(1,:)*curP1Grad(2,[2 1])';
       curP1Grad(2,:)*curP1Grad(3,[2 1])';
       curP1Grad(1,:)*curP1Grad(3,[2 1])'];       
D2U_h(2,1) = D2U_h(1,2);

D2U_h = reshape( (D2U_h(:)*ones(1,length(x)))',[length(x) 2 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigma_h = getSigma_h(x,y,curElem,curKappa,lvl,p)

gradU_h = p.statics.gradU_h;

curGradU_h =  gradU_h(x,y,curElem,lvl,p);
curGradU_h = reshape (curGradU_h',2,1,[]);
sigma_h = matMul(curKappa,curGradU_h);
sigma_h = reshape(sigma_h,2,[])';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basisU = getDisplacementBasis(x,y,curElem,lvl,p)

basisU = P2Basis(x,y,curElem,lvl,p)*p.statics.basisCoefficients';

function basisU = P2Basis(x,y,curElem,lvl,p)

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

basisU = [b1;
          b2;
          b3;
          b1.*b2;
          b2.*b3;
          b1.*b3]';

