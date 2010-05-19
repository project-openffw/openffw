function p = ODP2postProc(p)
% author: David Guenther 
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
% load computed solution
x = p.level(end).x;

% load enumerated data
dofU4e = p.level(end).enum.dofU4e;
nrNodes = p.level(end).nrNodes;
nrEdges = p.level(end).nrEdges;
nrElems = p.level(end).nrElems;

%% post-processing of calculated data
u = x(1:nrNodes+nrEdges);
u4e = zeros(nrElems,6);

for curElem = 1:nrElems
	localDoF = dofU4e(curElem,:);    
    curU = u(localDoF);
    u4e(curElem,:) = curU;
end

%% OUTPUT 
p.level(end).u = u;
p.level(end).u4e = u4e;
