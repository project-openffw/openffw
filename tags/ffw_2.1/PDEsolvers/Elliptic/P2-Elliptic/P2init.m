function p = P2init(p)
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

% set coefficients matrix
p.statics.basisCoefficients = ...
    [ 1 0 0 -2  0 -2 
      0 1 0 -2 -2  0
      0 0 1  0 -2 -2
      0 0 0  4  0  0
      0 0 0  0  4  0
      0 0 0  0  0  4 ] ;
  
% integration parameters
%  -> up to which polynomial degree shall integration be exact?
p.params.integrationDegrees.createLinSys.Stima = 1;
p.params.integrationDegrees.createLinSys.Dama = 2;
p.params.integrationDegrees.createLinSys.Mama = 2;
p.params.integrationDegrees.createLinSys.Rhs = 2;
p.params.integrationDegrees.createLinSys.Neumann = 2;
p.params.integrationDegrees.estimate.jumpTerm = 2;
p.params.integrationDegrees.estimate.volumeTerm = 4;
p.params.integrationDegrees.estimate.oscTerm = 2;
return


%% Basis Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = getGradBasis(pts,curElem,lvl,p)
nrPts = size(pts,1);
P1grad4e = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;

curP1Grad = P1grad4e(:,:,curElem);
curBasisU = getP2Basis(pts,curElem,lvl,p);

val  = zeros(6,2,nrPts);
for i = 1:nrPts
    curGradP2 = [ curP1Grad ;
        curP1Grad(1,:)*curBasisU(i,2) + curBasisU(i,1)*curP1Grad(2,:) ;
        curP1Grad(2,:)*curBasisU(i,3) + curBasisU(i,2)*curP1Grad(3,:) ;
        curP1Grad(1,:)*curBasisU(i,3) + curBasisU(i,1)*curP1Grad(3,:) ] ;
    val(:,:,i) = C * curGradP2;
end


function val = getD2Basis(pts,curElem,lvl,p)
curP1Grad = p.level(lvl).enum.P1grad4e(:,:,curElem);
C = p.statics.basisCoefficients;

val = zeros(3,6);
val(1,:) = C * [0;0;0;
       2*curP1Grad(1,1)*curP1Grad(2,1);
       2*curP1Grad(2,1)*curP1Grad(3,1);
       2*curP1Grad(1,1)*curP1Grad(3,1)];
val(2,:) = C * [0;0;0;
       2*curP1Grad(1,2)*curP1Grad(2,2);
       2*curP1Grad(2,2)*curP1Grad(3,2);
       2*curP1Grad(1,2)*curP1Grad(3,2)];
val(3,:) = C * [0;0;0;
       curP1Grad(1,:)*curP1Grad(2,[2 1])';
       curP1Grad(2,:)*curP1Grad(3,[2 1])';
       curP1Grad(1,:)*curP1Grad(3,[2 1])'];       
val = reshape( val(:)*ones(1,length(x)),[3 6 length(x)]);


function val = getBasis(pts,curElem,lvl,p)
val = getP2Basis(pts,curElem,lvl,p)*p.statics.basisCoefficients';


function val = getP2Basis(pts,curElem,lvl,p)
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

val = [b1;
          b2;
          b3;
          b1.*b2;
          b2.*b3;
          b1.*b3]';
      
      
%% Basis Functions - vectorized %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = getGradBasisVectorised(pts,pts_ref,parts,lvl,p)
nrPts = size(pts,1);
P1grad4e = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;
basisU = getP2BasisVectorised(pts,pts_ref,parts,lvl,p);

gradP2  = zeros(6,2,nrPts);
gradP2(1:3,:,:) = P1grad4e(:,:,parts);
P1grad4e = permute(P1grad4e,[ 3 1 2 ]);
gradP2(4,:,:) = ([P1grad4e(parts,1,1).*basisU(:,2),P1grad4e(parts,1,2).*basisU(:,2)]...
              + [basisU(:,1).*P1grad4e(parts,2,1),basisU(:,1).*P1grad4e(parts,2,2)])';
gradP2(5,:,:) = ([P1grad4e(parts,2,1).*basisU(:,3),P1grad4e(parts,2,2).*basisU(:,3)]...
              + [basisU(:,2).*P1grad4e(parts,3,1),basisU(:,2).*P1grad4e(parts,3,2)])';
gradP2(6,:,:) = ([P1grad4e(parts,1,1).*basisU(:,3),P1grad4e(parts,1,2).*basisU(:,3)]...
              + [basisU(:,1).*P1grad4e(parts,3,1),basisU(:,1).*P1grad4e(parts,3,2)])';

C = reshape(C(:)*ones(1,nrPts),[size(C,1) size(C,2) nrPts]);

val = matMul(C,gradP2);


function val = getD2BasisVectorised(pts,pts_ref,parts,lvl,p)
nrPts = size(pts,1);
P1grad4e = p.level(lvl).enum.P1grad4e;
C = p.statics.basisCoefficients;
C = reshape(C(:)*ones(1,nrPts),[size(C,1) size(C,2) nrPts]);
P1grad4e = permute(P1grad4e,[ 3 1 2 ]);

val = zeros(3,6,nrPts);

dummy = [zeros(nrPts,3),...
       2*P1grad4e(parts,1,1).*P1grad4e(parts,2,1),...
       2*P1grad4e(parts,2,1).*P1grad4e(parts,3,1),...
       2*P1grad4e(parts,1,1).*P1grad4e(parts,3,1)]';
dummy = reshape(dummy,[6 1 nrPts]);
val(1,:,:) = matMul(C,dummy);

dummy = [zeros(nrPts,3),...
       2*P1grad4e(parts,1,2).*P1grad4e(parts,2,2),...
       2*P1grad4e(parts,2,2).*P1grad4e(parts,3,2),...
       2*P1grad4e(parts,1,2).*P1grad4e(parts,3,2)]';
dummy = reshape(dummy,[6 1 nrPts]);
val(2,:,:) = matMul(C,dummy);

dummy = [zeros(nrPts,3),...
       P1grad4e(parts,1,1).*P1grad4e(parts,2,2)+P1grad4e(parts,1,2).*P1grad4e(parts,2,1),...
       P1grad4e(parts,2,1).*P1grad4e(parts,3,2)+P1grad4e(parts,2,2).*P1grad4e(parts,3,1),...
       P1grad4e(parts,1,1).*P1grad4e(parts,3,2)+P1grad4e(parts,1,2).*P1grad4e(parts,3,1)]';
dummy = reshape(dummy,[6 1 nrPts]);
val(3,:,:) = matMul(C,dummy);


function val = getBasisVectorised(pts,pts_ref,parts,lvl,p)
val = getP2BasisVectorised(pts,pts_ref,parts,lvl,p)*p.statics.basisCoefficients';


function val = getP2BasisVectorised(pts,pts_ref,parts,lvl,p)
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

val = [b1,...
       b2,...
       b3,...
       b1.*b2,...
       b2.*b3,...
       b1.*b3];