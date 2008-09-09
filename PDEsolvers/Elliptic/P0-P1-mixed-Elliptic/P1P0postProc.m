function p = P1P0postProc(p)
% computes the post-processing datas

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


%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;

% load discrete solution 
x = p.level(end).x;

% load enumerated data
nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = x(2*nrNodes+1:end)';

pU = x(1:nrNodes)';
pV = x(nrNodes+1:2*nrNodes)';

pU = pU(n4e);
pV = pV(n4e);

grad = zeros(3,2,nrElems);
grad(:,1,:) = pU';
grad(:,2,:) = pV';

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).u = u;
p.level(end).grad = grad;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
