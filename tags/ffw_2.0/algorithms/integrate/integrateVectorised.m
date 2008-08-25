function val = integrateVectorised( parts, curLvl,  degree, integrand, p )
% Integrates function handle integrand over given domain 1D or 2D.
% Using an outer loop over the Gauss points.
% 'integrand' must have the form integrand(pts,pts_ref,parts,lvl,p) and the
% return value must have the dimension [ n m nrElems ].

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


%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e4n = p.level(curLvl).enum.e4n;
ed4n = p.level(curLvl).enum.ed4n;
nrParts = size(parts,1);

%% Number of Gauss Points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrGaussPoints = ceil((degree+1)/2);

%% 1D or 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemType = size(parts,2);

if elemType == 2
% 1D, i.e. edge
    [refGaussPoints,refGaussWeights] = getGaussPoints(nrGaussPoints);
    refGaussPoints = [refGaussPoints,zeros(nrGaussPoints,1)];
    indices = rowaddr(ed4n,parts(:,1),parts(:,2));
    weightFactor = p.statics.getDetDPhiEd(indices,refGaussPoints,curLvl,p);
    
elseif elemType == 3 
% 2D triangle    
    [refGaussPoints,refGaussWeights] = getConProdGaussPoints(nrGaussPoints);
    indices = rowaddr(e4n,parts(:,1),parts(:,2));
    weightFactor = p.statics.getDetDPhiE(indices,refGaussPoints,curLvl,p);
    nrGaussPoints = nrGaussPoints^2;

elseif elemType == 4 
% 2D square
    [refGaussPoints,refGaussWeights] = getCompositeGaussPoints(nrGaussPoints);
    indices = rowaddr(e4n,parts(:,1),parts(:,2));
    weightFactor = p.statics.getDetDPhiE(indices,refGaussPoints,curLvl,p);
    nrGaussPoints = nrGaussPoints^2;
    
else
    error('Could not integrate because no Domain is specified.');
end

%% Integrate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for curGaussPt = 1:nrGaussPoints
    % transformation
    curCoords4e = p.statics.ref2arbitrary(parts,refGaussPoints(curGaussPt,:),curLvl,p);
    % evaluate function handle eval must be [n m nrParts]
    integrandVal = integrand(curCoords4e, refGaussPoints(curGaussPt,:), ...
                             indices, curLvl, p);
    if(curGaussPt == 1)
        dim = size(integrandVal);
        val = zeros(dim(1),dim(2),nrParts);
    end
    % weighted sum
    curWeights4e = weightFactor(:,curGaussPt)*refGaussWeights(curGaussPt);    
    curWeights4e = reshape( (curWeights4e*ones(1,dim(1)*dim(2)))',[dim(1),dim(2),nrParts]);
    val = val + curWeights4e.*integrandVal;
end

%% Result %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val = permute(val,[ 3 1 2 ]);