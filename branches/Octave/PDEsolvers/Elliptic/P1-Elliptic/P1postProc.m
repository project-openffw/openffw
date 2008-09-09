function p = P1postProc(p)
% computes the post-processing datas
% for a conforming P1-FE method.

% Copyright 2007 Jan Reininghaus, David Guenther, Joscha Gedicke
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

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = p.level(end).x;

% load enumerated data
grad4e = p.level(end).enum.grad4e;
nrNodes = p.level(end).nrNodes;
dofU4e = p.level(end).enum.dofU4e;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = x(1:nrNodes);

u4e = u(dofU4e);

grad_x = permute(grad4e(:,1,:),[3 1 2]);
grad_y = permute(grad4e(:,2,:),[3 1 2]);
gradU4e = [sum(u4e.*grad_x,2),sum(u4e.*grad_y,2)];

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).u = u;
p.level(end).u4e = u4e;
p.level(end).gradU4e = gradU4e;

