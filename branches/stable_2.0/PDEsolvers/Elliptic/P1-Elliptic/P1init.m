function p = P1init(p)
% makes available all necessary initial data

% Copyright 2007 Jan Reininghaus, David Guenther, 
%                Joscha Gedicke, Andreas Byfut
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

% integration parameters
%  -> up to which polynomial degree shall integration be exact?
p.params.integrationDegrees.createLinSys.Stima = 1;
p.params.integrationDegrees.createLinSys.Dama = 1;
p.params.integrationDegrees.createLinSys.Mama = 1;
p.params.integrationDegrees.createLinSys.Rhs = 1;
p.params.integrationDegrees.createLinSys.Neumann = 1;
p.params.integrationDegrees.estimate.jumpTerm = 1;
p.params.integrationDegrees.estimate.volumeTerm = 2;
p.params.integrationDegrees.estimate.oscTerm = 1;
return


%% Basis Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = getBasis(pts,curElem,lvl,p)
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
val = [b1;b2;b3]';


function val = getGradBasis(pts,curElem,lvl,p)
nrPts = size(pts,1);
val = p.level(lvl).enum.grad4e(:,:,curElem);
val = reshape(val(:)*ones(1,nrPts),[3 2 nrPts]);


function val = getD2Basis(pts,curElem,lvl,p)
nrPts = size(pts,1);
val = zeros(3,3,nrPts);


%% Basis Functions - vectorized %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = getBasisVectorised(pts,pts_ref,parts,lvl,p)
x = pts(:,1);
y = pts(:,2);
x_ref = pts_ref(:,1);
y_ref = pts_ref(:,2);

if y_ref~= 0
val = [(1-x_ref-y_ref)*ones(length(x),1), ...
           x_ref*ones(length(x),1), ...
           y_ref*ones(length(x),1) ];
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
val = [b1,b2,b3];    
end


function val = getGradBasisVectorised(pts,pts_ref,parts,lvl,p)
val = p.level(lvl).enum.grad4e(:,:,parts);


function val = getD2BasisVectorised(pts,pts_ref,parts,lvl,p)
nrPts = size(pts,1); 
val = zeros(3,3,nrPts);