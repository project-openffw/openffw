% Copyright 2007 David Guenther
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
%
%
%

clear all
syms x y real


u = sin(x^3)*cos(y^3);
% u = r.^(2/3).*sin(2/3*theta) + eps2_osc*cos(theta)*sin(eps_osc+k*r_osc)...
%         .*((r_osc/eps_osc).^4-2*(r_osc/eps_osc).^2+1);

% Volume force
matrix = [2 0; 0 2];
      
gradU = matrix*[diff(u,x); diff(u,y)]
KappaGradU = gradU;
DiffKappaGradU = [diff(KappaGradU(1),x) diff(KappaGradU(2),x);
                  diff(KappaGradU(1),y) diff(KappaGradU(2),y)]
DiffKappaGradU = matrix*DiffKappaGradU;  
minusDivKappaGradU = -simple(DiffKappaGradU(1,1) + DiffKappaGradU(2,2));



f = minusDivKappaGradU;
f = Matlab4Maple(f)
