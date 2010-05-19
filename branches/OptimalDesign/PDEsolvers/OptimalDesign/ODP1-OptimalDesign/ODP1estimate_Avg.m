function p = ODP1estimate_Avg(p)
%author: David Guenther
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

%% INPUT
n4e = p.level(end).geom.n4e;
lvl = size(p.level,2);

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);
%% compute the average error 
eta4T = integrate(n4e,lvl,degree,@integrand,p);

%% OUTPUT
p.level(end).etaT = sqrt(eta4T);
p.level(end).estimatedError = sqrt(sum(eta4T));

%% supply the integrand ||p_h - Ap_h||_L^2
function val = integrand(points,curElem,lvl,p)

Ap_h = p.statics.Ap_h;
sigma_h = p.statics.sigma_h;

Aph = Ap_h(points,curElem,lvl,p);
sigmah = sigma_h(points,curElem,lvl,p);

val = sum((Aph - sigmah).*(Aph - sigmah),2);

val = reshape(val,[1 1 length(points(:,1))]);
    
