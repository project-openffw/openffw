function drawDisplacementBasis(j,res,p)
% draws the basis functions of an finite element

% Copyright 2007 Jan Reininghaus, David Guenther, 
%                Andreas Byfut, Joscha Gedicke
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


%% input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = loadField('p.statics','basis',p,[]);
if isempty(u)
    u = p.statics.basisU;
end
nrBasisFunc = size(p.level(end).enum.dofU4e,2);

%% draw routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XI = 0:1/res:1;
YI = 0:1/res:1;

UI = zeros(res);

for indX = 1:length(XI)
	for indY = 1:length(YI)
        x = XI(indX);
        y = YI(indY);
        if(x+y > 1)
            val = NaN*ones(1,nrBasisFunc);
        else
            val  = u([x,y],1,1,p);
        end
        UI(indX,indY) = val(j);
	end
end

surf(XI,YI,UI','EdgeColor','none');
