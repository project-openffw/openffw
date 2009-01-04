function p = RT0P0init(p)
% makes available all necessary initial data,
% handels the discrete displacement,flux, 

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
p.statics.basisU = @getDisplacementBasis;
p.statics.u_h = @getU_h;
p.statics.sigma_h = @getSigma_h;

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
function u_h = getU_h(pts,curElem,lvl,p)
nrPts = size(pts,1);
u = p.level(lvl).u4e;
curU = u(curElem);
u_h = repmat(curU,1,nrPts);


function sigma_h = getSigma_h(pts,curElem,curKappa,lvl,p)
nrPts = size(pts,1);
basisSigma = getStressBasis(pts,curElem,lvl,p);
grad4e = p.level(lvl).grad4e;
sigma_h = zeros(nrPts,2);
for j = 1:nrPts
    sigma_h(j,:) = curKappa(:,:,j)*(grad4e(curElem,:)*basisSigma(:,:,j))';
end


function basisSigma = getStressBasis(pts,curElem,lvl,p)
nrPts = size(pts,1);
x = pts(:,1);
y = pts(:,2);
c4n = p.level(lvl).geom.c4n;
n4e = p.level(lvl).geom.n4e;
ed4e = p.level(lvl).enum.ed4e;
e4ed = p.level(lvl).enum.e4ed;
area4e = p.level(lvl).enum.area4e;
length4ed = p.level(lvl).enum.length4ed;

edges = ed4e(curElem,:);
coords = c4n(n4e(curElem,:),:);
area = area4e(curElem);

P1 = coords(1,:)';
P2 = coords(2,:)';
P3 = coords(3,:)';

signum = ones(3,1);
I = find(e4ed(edges,2) == curElem);
signum(I) = -1;

basisSigma = zeros(3,2,nrPts);

for j = 1:nrPts
    curX = x(j);
    curY = y(j);
    
    b1 = signum(1)*length4ed(edges(1))/(2*area)*([curX;curY] - P3);	
    b2 = signum(2)*length4ed(edges(2))/(2*area)*([curX;curY] - P1);	
    b3 = signum(3)*length4ed(edges(3))/(2*area)*([curX;curY] - P2);

    basisSigma(:,:,j) = [b1';b2';b3'];
end


function basisU = getDisplacementBasis(pts,curElem,lvl,p)
nrPts = size(pts,1);
basisU = ones(nrPts,1);