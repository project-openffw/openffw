function val = integrateVectorised( parts, curLvl,  degree, integrand, p )
% Integrates function handle integrand over given domain 1D or 2D.
% Using an outer loop over the Gauss points.
% 'integrand' must have the form integrand(x,y,x_ref,y_ref,parts,lvl,p) and the
% return value must have the dimension [ n m nrElems ].

% Copyright 2007 Joscha Gedicke
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


%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
area4e = p.level(curLvl).enum.area4e;
length4ed = p.level(curLvl).enum.length4ed;
e4n = p.level(curLvl).enum.e4n;
ed4n = p.level(curLvl).enum.ed4n;
nrParts = size(parts,1);

%% Number of Gauss Points%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrGaussPoints = ceil((degree+1)/2);

%% 1D or 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(parts,2) == 2 
    % 1D
    [refGaussPoints,refGaussWeights] = getGaussPoints(nrGaussPoints);
    refGaussPoints = [refGaussPoints,zeros(nrGaussPoints,1)];
    indices = rowaddr(ed4n,parts(:,1),parts(:,2));
    weightFactor = length4ed(indices);
    
elseif size(parts,2) == 3 
    % 2D triangle
    [refGaussPoints,refGaussWeights] = getConProdGaussPoints(nrGaussPoints);
    indices = rowaddr(e4n,parts(:,1),parts(:,2));
    weightFactor = 2*area4e(indices);  
    nrGaussPoints = nrGaussPoints^2;

elseif size(parts,2) == 4 
    % 2D square
    [refGaussPoints,refGaussWeights] = getCompositeGaussPoints(nrGaussPoints);
    indices = rowaddr(e4n,parts(:,1),parts(:,2));
    weightFactor = area4e(indices);  
    nrGaussPoints = nrGaussPoints^2;
    
else
    error('Could not integrate because no Domain is specified.');
end

%% Integrate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for curPoint = 1:nrGaussPoints
    %transformation
    curCoords4e = ref2arbitrary(parts,refGaussPoints(curPoint,:),curLvl,p);
    % evaluate function handle eval must be [n m nrParts]
    integrandVal = integrand(curCoords4e(:,1), curCoords4e(:,2),...
                       refGaussPoints(curPoint,1), refGaussPoints(curPoint,2), ...
                       indices, curLvl, p);
    if(curPoint == 1)
        dim = size(integrandVal);
        val = zeros(dim(1),dim(2),nrParts);
    end
    % weighted sum
    curWeights4e = weightFactor*refGaussWeights(curPoint);
    curWeights4e = reshape( (curWeights4e*ones(1,dim(1)*dim(2)))',[dim(1),dim(2),nrParts]);
    val = val + curWeights4e.*integrandVal;
end

%% Result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val = permute(val,[ 3 1 2 ]);


%% ref2arbitrary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = ref2arbitrary(parts,curPoint,curLvl,p)

c4n = p.level(curLvl).geom.c4n;

%1D or 2D
if size(parts,2) == 2
    %1D ref = conv(0,1)
    c1 = c4n(parts(:,1),:);
    c2 = c4n(parts(:,2),:);
    
    val = c1 + curPoint(1)*(c2-c1);   

elseif size(parts,2) == 3 
    % 2D ref = conv{ (0,0),(1,0),(0,1)}
    c1 = c4n(parts(:,1),:);
    c2 = c4n(parts(:,2),:);
    c3 = c4n(parts(:,3),:);
    
    val = c1 + curPoint(1)*(c2-c1) + ...
               curPoint(2)*(c3-c1);

elseif size(parts,2) == 4
    % 2D ref = conv{(0,0),(1,0),(0,1),(1,1)}
    c1 = c4n(parts(:,1),:);
    c2 = c4n(parts(:,2),:);
    c3 = c4n(parts(:,3),:);
    c4 = c4n(parts(:,4),:);
    
    val = c1 + curPoint(1)*(c2-c1) + ...
               curPoint(2)*(c4-c1 + curPoint(1)*(c3-c4-c2+c1));
end
