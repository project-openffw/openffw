function p = Elasticity_Template(p)

% Copyright 2007 Jan Reininghaus
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


%% Problem definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDE definition
p.problem.geom = 'Square';
p.problem.f = @f;
p.problem.g = @g;
p.problem.u_D = @u_D;

% set Lamé parameters
p.PDE.E = 2900; 
p.PDE.nu = 0.3;
return

%% Volume force %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = f(pts,p)
z = ones(size(pts,1),2);

%% Dirichlet boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = u_D(pts,p)
z = zeros(size(pts,1),2);

%% Neumann boundary values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = g(pts,n,p)
z = zeros(size(pts,1),2);