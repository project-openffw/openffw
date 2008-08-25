function p = elasticity_setLame(p,nu,E)
% set lame parameters

% Copyright 2007 Jan Reininghaus, David Guenther
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


if nargin > 1
    p.PDE.nu = nu;
elseif nargin > 2
    p.PDE.nu = nu;
    p.PDE.E = E;
end

p.PDE.mu = p.PDE.E/(2*(1+p.PDE.nu));
p.PDE.lambda = p.PDE.E * p.PDE.nu /( (1+p.PDE.nu)*(1-2*p.PDE.nu) );

return




