function p = transformation_affin(p)
% This function initializes the transformation from the reference element
% to the local elements with a biaffin mapping

% Copyright 2007 Andreas Byfut, Joscha Gedicke
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


%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.statics.ref2arbitrary = @ref2arbitrary;
p.statics.getDetDPhiE = @getDetDPhiE;
p.statics.getDetDPhiEd = @getDetDPhiEd;
return


%% transformation functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = ref2arbitrary(n4parts,refPts,curLvl,p)
c4n = p.level(curLvl).geom.c4n;
nrRefPts = size(refPts,1);
nrParts = size(n4parts,1);
dimParts = size(n4parts,2);

n4parts = repmat(n4parts,[1 1 nrRefPts]);
n4parts = reshape(permute(n4parts,[3 1 2]),[nrParts*nrRefPts dimParts]);
refPts = repmat(refPts,[nrParts 1]);

%1D or 2D
if size(n4parts,2) == 2
    %1D ref = conv(0,1)
    c1 = c4n(n4parts(:,1),:);
    c2 = c4n(n4parts(:,2),:);    
    val = c1 + [refPts(:,1),refPts(:,1)].*(c2-c1);   

elseif size(n4parts,2) == 3 
    % 2D ref = conv{ (0,0),(1,0),(0,1)}
    c1 = c4n(n4parts(:,1),:);
    c2 = c4n(n4parts(:,2),:);
    c3 = c4n(n4parts(:,3),:);
    val = c1 + [refPts(:,1),refPts(:,1)].*(c2-c1) + ...
               [refPts(:,2),refPts(:,2)].*(c3-c1);

elseif size(n4parts,2) == 4
    % 2D ref = conv{(0,0),(1,0),(0,1),(1,1)}
    c1 = c4n(n4parts(:,1),:);
    c2 = c4n(n4parts(:,2),:);
    c3 = c4n(n4parts(:,3),:);
    c4 = c4n(n4parts(:,4),:);
    val = c1 + [refPts(:,1),refPts(:,1)].*(c2-c1) + ...
               [refPts(:,2),refPts(:,2)].*(c4-c1) + ...
               [refPts(:,1),refPts(:,1)].*[refPts(:,2),refPts(:,2)].* ...
                    (c3-c4-c2+c1);
end



function val = getDetDPhiE(parts,refPts,curLvl,p)
% parts = indicies for parts
nrRefPts = size(refPts,1);
val = 2*p.level(curLvl).enum.area4e(parts);
val = repmat(val,[1 nrRefPts]);


function val = getDetDPhiEd(parts,refPts,curLvl,p)
% parts = indicies for parts
nrRefPts = size(refPts,1);
val = p.level(curLvl).enum.length4ed(parts);
val = repmat(val,[1 nrRefPts]);