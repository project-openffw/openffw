function p = P3init(p)
% makes available all necessary initial data

% Copyright 2007 Joscha Gedicke, Andreas Byfut
%
% This file is part of FFW.
%
% FFW is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% FFW is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set function handles
p.statics.basis = @getBasis;
p.statics.gradBasis = @getGradBasis;
p.statics.d2Basis = @getD2Basis;

% set function handles vectorised
p.statics.basisVectorised = @getBasisVectorised;
p.statics.gradBasisVectorised = @getGradBasisVectorised;
p.statics.d2BasisVectorised = @getD2BasisVectorised;

p.statics.basisCoefficients = P3getBasisCoefficients();

% integration parameters 
%  -> up to which polynomial degree shall integration be exact?
p.params.integrationDegrees.createLinSys.Stima = 4;
p.params.integrationDegrees.createLinSys.Dama = 5;
p.params.integrationDegrees.createLinSys.Mama = 6;
p.params.integrationDegrees.createLinSys.Rhs = 6;
p.params.integrationDegrees.createLinSys.Neumann = 6;
p.params.integrationDegrees.estimate.jumpTerm = 4;
p.params.integrationDegrees.estimate.volumeTerm = 6;
p.params.integrationDegrees.estimate.oscTerm = 6;
return


%% Basis Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = getGradBasis(pts,curElem,lvl,p)
nrPts = size(pts,1);
P1grad4e = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;

curP1Grad = P1grad4e(:,:,curElem);
curBasisU = getP3Basis(pts,curElem,lvl,p);

val = zeros(10,2,nrPts);

for i = 1 : nrPts
    u = curBasisU(i,:)';
    curGradP3 = [ curP1Grad ;
        curP1Grad(1,:)*u(2) + u(1)*curP1Grad(2,:) ;
        curP1Grad(2,:)*u(3) + u(2)*curP1Grad(3,:) ;
        curP1Grad(3,:)*u(1) + u(3)*curP1Grad(1,:) ;
        [u(2)*u(1),u(1)*u(1),u(1)*u(2)]*curP1Grad([1 2 1],:) - [u(2)*u(2),u(1)*u(2),u(1)*u(2)]*curP1Grad([1 2 2],:);
        [u(3)*u(2),u(2)*u(2),u(2)*u(3)]*curP1Grad([2 3 2],:) - [u(3)*u(3),u(2)*u(3),u(2)*u(3)]*curP1Grad([2 3 3],:);
        [u(1)*u(3),u(3)*u(3),u(3)*u(1)]*curP1Grad([3 1 3],:) - [u(1)*u(1),u(3)*u(1),u(3)*u(1)]*curP1Grad([3 1 1],:);
        [u(2)*u(3),u(1)*u(3),u(1)*u(2)]*curP1Grad ] ;
    val(:,:,i) = C * curGradP3;
end


function val = getD2Basis(pts,curElem,lvl,p)
nrPts = size(pts,1);
P1grad4e = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;
curP1Grad = P1grad4e(:,:,curElem);
% curBasis = getP3Basis(pts,curElem,lvl,p);
val = zeros(3,10,nrPts);

for i = 1 : nrPts  
D2U_h(1,:,:) = C * [0;0;0;
       2*curP1Grad(1,1)*curP1Grad(2,1);
       2*curP1Grad(2,1)*curP1Grad(3,1);
       2*curP1Grad(1,1)*curP1Grad(3,1);
       d2Help(curP1Grad,u,1,2,1,1,1) - d2Help(curP1Grad,u,1,2,2,1,1);
       d2Help(curP1Grad,u,2,3,2,1,1) - d2Help(curP1Grad,u,2,3,3,1,1);
       d2Help(curP1Grad,u,3,1,3,1,1) - d2Help(curP1Grad,u,3,1,1,1,1);
       d2Help(curP1Grad,u,1,2,3,1,1)];
D2U_h(2,:,:) = C * [0;0;0;
       2*curP1Grad(1,2)*curP1Grad(2,2);
       2*curP1Grad(2,2)*curP1Grad(3,2);
       2*curP1Grad(1,2)*curP1Grad(3,2);
       d2Help(curP1Grad,u,1,2,1,2,2) - d2Help(curP1Grad,u,1,2,2,2,2);
       d2Help(curP1Grad,u,2,3,2,2,2) - d2Help(curP1Grad,u,2,3,3,2,2);
       d2Help(curP1Grad,u,3,1,3,2,2) - d2Help(curP1Grad,u,3,1,1,2,2);
       d2Help(curP1Grad,u,1,2,3,2,2)];
D2U_h(3,:,:) = C * [0;0;0;
       curP1Grad(1,:)*curP1Grad(2,[2 1])';
       curP1Grad(2,:)*curP1Grad(3,[2 1])';
       curP1Grad(1,:)*curP1Grad(3,[2 1])';
       d2Help(curP1Grad,u,1,2,1,1,2) - d2Help(curP1Grad,u,1,2,2,1,2);
       d2Help(curP1Grad,u,2,3,2,1,2) - d2Help(curP1Grad,u,2,3,3,1,2);
       d2Help(curP1Grad,u,3,1,3,1,2) - d2Help(curP1Grad,u,3,1,1,1,2);
       d2Help(curP1Grad,u,1,2,3,1,2)];
end


function val = d2Help(grad,u,i,j,k,d1,d2)
val = [ ...
    grad(j,d2)*u(k) + u(j)*grad(k,d2) , ...
    grad(k,d2)*u(i) + u(k)*grad(i,d2) , ...
    grad(i,d2)*u(j) + u(i)*grad(j,d2) ] * ...
      grad([i j k],d1);


function val = getBasis(pts,curElem,lvl,p)
val = getP3Basis(pts,curElem,lvl,p)*p.statics.basisCoefficients';


function val = getP3Basis(pts,curElem,lvl,p)
x = pts(:,1);
y = pts(:,2);
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

val =    [b1;
          b2;
          b3;
          b1.*b2;
          b2.*b3;
          b3.*b1;
          b1.*b2.*(b1-b2);
          b2.*b3.*(b2-b3);
          b3.*b1.*(b3-b1);
          b1.*b2.*b3 ]';
      
      
%% Basis Functions - vectorized %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = getGradBasisVectorised(pts,pts_ref,parts,lvl,p)
nrPts = size(pts,1);
P1grad4e = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;
basis = getP3BasisVectorised(pts,pts_ref,parts,lvl,p);

val  = zeros(10,2,nrPts);

val(1:3,:,:) = P1grad4e(:,:,parts);

P1grad4e = permute(P1grad4e,[ 3 1 2 ]);

val(4,:,:) = ([P1grad4e(parts,1,1).*basis(:,2),P1grad4e(parts,1,2).*basis(:,2)]...
              + [basis(:,1).*P1grad4e(parts,2,1),basis(:,1).*P1grad4e(parts,2,2)])';
val(5,:,:) = ([P1grad4e(parts,2,1).*basis(:,3),P1grad4e(parts,2,2).*basis(:,3)]...
              + [basis(:,2).*P1grad4e(parts,3,1),basis(:,2).*P1grad4e(parts,3,2)])';
val(6,:,:) = ([P1grad4e(parts,1,1).*basis(:,3),P1grad4e(parts,1,2).*basis(:,3)]...
              + [basis(:,1).*P1grad4e(parts,3,1),basis(:,1).*P1grad4e(parts,3,2)])';
          
P1grad4e = permute(P1grad4e,[ 2 3 1 ]);
P1grad4e = P1grad4e(:,:,parts);

val(7,:,:) = gradP3Help(P1grad4e,basis,1,2,1)-gradP3Help(P1grad4e,basis,1,2,2);
val(8,:,:) = gradP3Help(P1grad4e,basis,2,3,2)-gradP3Help(P1grad4e,basis,2,3,3);
val(9,:,:) = gradP3Help(P1grad4e,basis,3,1,3)-gradP3Help(P1grad4e,basis,3,1,1);
val(10,:,:) = gradP3Help(P1grad4e,basis,1,2,3);

C = reshape(C(:)*ones(1,nrPts),[size(C,1) size(C,2) nrPts]);

val = matMul(C,val);

function val = gradP3Help(grad,u,i,j,k)
    val = [u(:,j).*u(:,k),u(:,i).*u(:,k),u(:,i).*u(:,j)];
    val = reshape(val',1,3,[]);
    val = matMul(val,grad([i j k],:,:));
    val = permute(val, [ 2 3 1 ]);
    

function val = getD2BasisVectorised(pts,pts_ref,parts,lvl,p)
nrPts = size(pts,1);
grad = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;
C = reshape(C(:)*ones(1,nrPts),[size(C,1) size(C,2) nrPts]);
P1grad4e = permute(grad,[ 3 1 2 ]);
curP1Grad= grad(:,:,parts);

u = getP3BasisVectorised(pts,pts_ref,parts,lvl,p);

val = zeros(3,10,nrPts);

dummy = [zeros(nrPts,3),...
       2*P1grad4e(parts,1,1).*P1grad4e(parts,2,1),...
       2*P1grad4e(parts,2,1).*P1grad4e(parts,3,1),...
       2*P1grad4e(parts,1,1).*P1grad4e(parts,3,1),...
       d2HelpV(curP1Grad,P1grad4e,u,1,2,1,1,1) - d2HelpV(curP1Grad,P1grad4e,u,1,2,2,1,1),...
       d2HelpV(curP1Grad,P1grad4e,u,2,3,2,1,1) - d2HelpV(curP1Grad,P1grad4e,u,2,3,3,1,1),...
       d2HelpV(curP1Grad,P1grad4e,u,3,1,3,1,1) - d2HelpV(curP1Grad,P1grad4e,u,3,1,1,1,1),...
       d2HelpV(curP1Grad,P1grad4e,u,1,2,3,1,1)]';
dummy = reshape(dummy,[10 1 nrPts]);
val(1,:,:) = matMul(C,dummy);

dummy = [zeros(nrPts,3),...
       2*P1grad4e(parts,1,2).*P1grad4e(parts,2,2),...
       2*P1grad4e(parts,2,2).*P1grad4e(parts,3,2),...
       2*P1grad4e(parts,1,2).*P1grad4e(parts,3,2),...
       d2HelpV(curP1Grad,P1grad4e,u,1,2,1,2,2) - d2HelpV(curP1Grad,P1grad4e,u,1,2,2,2,2),...
       d2HelpV(curP1Grad,P1grad4e,u,2,3,2,2,2) - d2HelpV(curP1Grad,P1grad4e,u,2,3,3,2,2),...
       d2HelpV(curP1Grad,P1grad4e,u,3,1,3,2,2) - d2HelpV(curP1Grad,P1grad4e,u,3,1,1,2,2),...
       d2HelpV(curP1Grad,P1grad4e,u,1,2,3,2,2)]';
dummy = reshape(dummy,[10 1 nrPts]);
val(2,:,:) = matMul(C,dummy);

dummy = [zeros(nrPts,3),...
       P1grad4e(parts,1,1).*P1grad4e(parts,2,2)+P1grad4e(parts,1,2).*P1grad4e(parts,2,1),...
       P1grad4e(parts,2,1).*P1grad4e(parts,3,2)+P1grad4e(parts,2,2).*P1grad4e(parts,3,1),...
       P1grad4e(parts,1,1).*P1grad4e(parts,3,2)+P1grad4e(parts,1,2).*P1grad4e(parts,3,1),...
       d2HelpV(curP1Grad,P1grad4e,u,1,2,1,1,2) - d2HelpV(curP1Grad,P1grad4e,u,1,2,2,1,2),...
       d2HelpV(curP1Grad,P1grad4e,u,2,3,2,1,2) - d2HelpV(curP1Grad,P1grad4e,u,2,3,3,1,2),...
       d2HelpV(curP1Grad,P1grad4e,u,3,1,3,1,2) - d2HelpV(curP1Grad,P1grad4e,u,3,1,1,1,2),...
       d2HelpV(curP1Grad,P1grad4e,u,1,2,3,1,2)]';
dummy = reshape(dummy,[10 1 nrPts]);
val(3,:,:) = matMul(C,dummy);


function val = d2HelpV(grad,grad2,u,i,j,k,d1,d2)
val = [ ...
    grad2(:,j,d2).*u(:,k) + u(:,j).*grad2(:,k,d2) , ...
    grad2(:,k,d2,:).*u(:,i) + u(:,k).*grad2(:,i,d2) , ...
    grad2(:,i,d2,:).*u(:,j) + u(:,i).*grad2(:,j,d2) ]; % * ...   grad([i j k],d1);
val = reshape(val',1,3,[]);
val = matMul(val,grad([i j k],d1,:));
val = permute(val, [ 3 2 1 ]);
  

function val = getBasisVectorised(pts,pts_ref,parts,lvl,p)
val = getP3BasisVectorised(pts,pts_ref,parts,lvl,p)*p.statics.basisCoefficients';


function val = getP3BasisVectorised(pts,pts_ref,parts,lvl,p)
nrPts = size(pts,1);
x = pts(:,1);
y = pts(:,2);
x_ref = pts_ref(:,1);
y_ref = pts_ref(:,2);

if y_ref~= 0
    b1 = (1-x_ref-y_ref)*ones(nrPts,1);
    b2 =          x_ref*ones(nrPts,1);
    b3 =           y_ref*ones(nrPts,1) ;

else
    n4e = p.level(lvl).geom.n4e;
    c4n = p.level(lvl).geom.c4n;
    area4e = p.level(lvl).enum.area4e(parts);

    P1 = c4n(n4e(parts,1),:);
    P2 = c4n(n4e(parts,2),:);
    P3 = c4n(n4e(parts,3),:);

    b1 = 0.5./area4e.*( (P2(:,2)-P3(:,2)).*x + (P3(:,1)-P2(:,1)).*y + P2(:,1).*P3(:,2)-P3(:,1).*P2(:,2) );
    b2 = 0.5./area4e.*( (P3(:,2)-P1(:,2)).*x + (P1(:,1)-P3(:,1)).*y + P3(:,1).*P1(:,2)-P1(:,1).*P3(:,2) );
    b3 = 0.5./area4e.*( (P1(:,2)-P2(:,2)).*x + (P2(:,1)-P1(:,1)).*y + P1(:,1).*P2(:,2)-P2(:,1).*P1(:,2) );
end

val = [b1, ...
          b2, ...
          b3, ...
          b1.*b2, ...
          b2.*b3, ...
          b1.*b3, ...
          b1.*b2.*(b1-b2), ...
          b2.*b3.*(b2-b3), ...
          b3.*b1.*(b3-b1), ...
          b1.*b2.*b3];