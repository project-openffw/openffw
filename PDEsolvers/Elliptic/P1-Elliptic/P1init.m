function p = P1init(p)
% makes available all necessary initial data

% Copyright 2007 Jan Reininghaus, David Guenther, Joscha Gedicke
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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = getBasis(x,y,curElem,lvl,p)

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

%%
function val = getGradBasis(x,y,curElem,lvl,p)

val = p.level(lvl).enum.grad4e(:,:,curElem);
val = reshape(val(:)*ones(1,length(x)),[3 2 length(x)]);

%%
function val = getD2Basis(x,y,curElem,lvl,p)

val = zeros(3,3,length(x));

%% Vectorised  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function val = getBasisVectorised(x,y,x_ref,y_ref,parts,lvl,p)

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

%%
function val = getGradBasisVectorised(x,y,x_ref,y_ref,parts,lvl,p)

val = p.level(lvl).enum.grad4e(:,:,parts);

%%
function val = getD2BasisVectorised(x,y,x_ref,y_ref,parts,lvl,p)

val = zeros(3,3,length(x));
