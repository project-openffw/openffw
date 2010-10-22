function p = P3run(p)

% author: Joscha Gedicke

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial enumeration
enumerate = p.statics.enumerate;
p = enumerate(p);

% Set initial values
nrNodes = p.level(1).nrNodes;
nrEdges = p.level(1).nrEdges;
nrElems = p.level(1).nrElems;
p.level(1).x = zeros(nrNodes+2*nrEdges+nrElems,1);

p.statics.basisU = @getDisplacementBasis;
p.statics.basisGradU = @getGradP3;
p.statics.u_h = @getU_h;
p.statics.gradU_h = @getGradU_h;
p.statics.D2U_h = @getD2U_h;
p.statics.sigma_h = @getSigma_h;

p.statics.basisCoefficients = P3getBasisCoefficients();

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

curGradP3 = getGradP3(x,y,curElem,lvl,p);

for i = 1 : length(x)
    gradU_h(i,:) = curU  * curGradP3(:,:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gradP3 = getGradP3(x,y,curElem,lvl,p)

P1grad4e = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;

curP1Grad = P1grad4e(:,:,curElem);
curBasisU = P3Basis(x,y,curElem,lvl,p);

gradP3 = zeros(10,2,length(x));

for i = 1 : length(x)
    u = curBasisU(i,:)';
    curGradP3 = [ curP1Grad ;
        curP1Grad(1,:)*u(2) + u(1)*curP1Grad(2,:) ;
        curP1Grad(2,:)*u(3) + u(2)*curP1Grad(3,:) ;
        curP1Grad(3,:)*u(1) + u(3)*curP1Grad(1,:) ;
        [u(2)*u(1),u(1)*u(1),u(1)*u(2)]*curP1Grad([1 2 1],:) - [u(2)*u(2),u(1)*u(2),u(1)*u(2)]*curP1Grad([1 2 2],:);
        [u(3)*u(2),u(2)*u(2),u(2)*u(3)]*curP1Grad([2 3 2],:) - [u(3)*u(3),u(2)*u(3),u(2)*u(3)]*curP1Grad([2 3 3],:);
        [u(1)*u(3),u(3)*u(3),u(3)*u(1)]*curP1Grad([3 1 3],:) - [u(1)*u(1),u(3)*u(1),u(3)*u(1)]*curP1Grad([3 1 1],:);
        [u(2)*u(3),u(1)*u(3),u(1)*u(2)]*curP1Grad ] ;
    gradP3(:,:,i) = C * curGradP3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D2U_h = getD2U_h(x,y,curElem,lvl,p)
% 
u = p.level(lvl).u4e;
P1grad4e = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;

curU = u(curElem,:);
curP1Grad = P1grad4e(:,:,curElem);
curBasisU = P3Basis(x,y,curElem,lvl,p);

D2U_h = zeros(length(x),2,2);

for i = 1 : length(x)
u = curBasisU(i,:);    
D2U_h(i,1,1) = curU * C * [0;0;0;
       2*curP1Grad(1,1)*curP1Grad(2,1);
       2*curP1Grad(2,1)*curP1Grad(3,1);
       2*curP1Grad(1,1)*curP1Grad(3,1);
       d2Help(curP1Grad,u,1,2,1,1,1) - d2Help(curP1Grad,u,1,2,2,1,1);
       d2Help(curP1Grad,u,2,3,2,1,1) - d2Help(curP1Grad,u,2,3,3,1,1);
       d2Help(curP1Grad,u,3,1,3,1,1) - d2Help(curP1Grad,u,3,1,1,1,1);
       d2Help(curP1Grad,u,1,2,3,1,1)];
D2U_h(i,2,2) = curU * C * [0;0;0;
       2*curP1Grad(1,2)*curP1Grad(2,2);
       2*curP1Grad(2,2)*curP1Grad(3,2);
       2*curP1Grad(1,2)*curP1Grad(3,2);
       d2Help(curP1Grad,u,1,2,1,2,2) - d2Help(curP1Grad,u,1,2,2,2,2);
       d2Help(curP1Grad,u,2,3,2,2,2) - d2Help(curP1Grad,u,2,3,3,2,2);
       d2Help(curP1Grad,u,3,1,3,2,2) - d2Help(curP1Grad,u,3,1,1,2,2);
       d2Help(curP1Grad,u,1,2,3,2,2)];
D2U_h(i,1,2) = curU * C * [0;0;0;
       curP1Grad(1,:)*curP1Grad(2,[2 1])';
       curP1Grad(2,:)*curP1Grad(3,[2 1])';
       curP1Grad(1,:)*curP1Grad(3,[2 1])';
       d2Help(curP1Grad,u,1,2,1,1,2) - d2Help(curP1Grad,u,1,2,2,1,2);
       d2Help(curP1Grad,u,2,3,2,1,2) - d2Help(curP1Grad,u,2,3,3,1,2);
       d2Help(curP1Grad,u,3,1,3,1,2) - d2Help(curP1Grad,u,3,1,1,1,2);
       d2Help(curP1Grad,u,1,2,3,1,2)];
D2U_h(i,2,1) = D2U_h(i,1,2);

end

function val = d2Help(grad,u,i,j,k,d1,d2)
val = [ ...
    grad(j,d2)*u(k) + u(j)*grad(k,d2) , ...
    grad(k,d2)*u(i) + u(k)*grad(i,d2) , ...
    grad(i,d2)*u(j) + u(i)*grad(j,d2) ] * ...
      grad([i j k],d1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigma_h = getSigma_h(x,y,curElem,curKappa,lvl,p)

gradU_h = p.statics.gradU_h;

curGradU_h =  gradU_h(x,y,curElem,lvl,p);
curGradU_h = reshape (curGradU_h',2,1,[]);
sigma_h = matMul(curKappa,curGradU_h);
sigma_h = reshape(sigma_h,2,[])';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function basisU = getDisplacementBasis(x,y,curElem,lvl,p)

basisU = P3Basis(x,y,curElem,lvl,p)*p.statics.basisCoefficients';

function basisU = P3Basis(x,y,curElem,lvl,p)

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
          b3.*b1;
          b1.*b2.*(b1-b2);
          b2.*b3.*(b2-b3);
          b3.*b1.*(b3-b1);
          b1.*b2.*b3 ]';
