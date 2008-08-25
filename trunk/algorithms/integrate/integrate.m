function val = integrate( parts, curLvl,  degree, integrand, p )
% Integrates function handle integrand over given domain 1D or 2D.
% Using an outer loop over the elements.
% 'integrand' must have the form integrand(pts,curPart,curLvl,p) and the
% return value must have the dimension [ n m nrGaussPoints ]

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

%% 1D or = 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemType = size(parts,2);
if elemType == 2 
    % 1D
    [refGaussPoints,refGaussWeights] = getGaussPoints(nrGaussPoints);
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
for curPart = 1:nrParts   
    % transformation
    curCoords = p.statics.ref2arbitrary(parts(curPart,:),refGaussPoints,curLvl,p);
    % evaluate function handle eval must be [n m nrGaussPoints]
    integrandVal = integrand(curCoords, indices(curPart), curLvl, p);
    if(curPart == 1)
        dim = size(integrandVal);
        val = zeros(nrParts,dim(1),dim(2));
    end
    % weighted sum
    curWeights = weightFactor(curPart,:)'.*refGaussWeights;
    curWeights = reshape( (curWeights*ones(1,dim(1)*dim(2)))',[dim(1),dim(2),nrGaussPoints]);
    val(curPart,:,:) =  sum(curWeights.*integrandVal,3);
    
end