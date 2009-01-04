function p = RT0P0postProc(p)
% computes the post-processing datas for 
% the mixed RT0-P0-FE method.

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


%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;

x = p.level(end).x;

% load enumerated data
ed4e = p.level(end).enum.ed4e;
e4ed = p.level(end).enum.e4ed;
length4ed = p.level(end).enum.length4ed;
area4e = p.level(end).enum.area4e;
nrElems = p.level(end).nrElems;
nrEdges = p.level(end).nrEdges;

%% PostProc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = x(nrEdges+1:end)';
gradP = x(1:nrEdges);

grad = zeros(3,2,nrElems);
grad4e = zeros(nrElems,3);

for curElem = 1:nrElems
	curEdges = ed4e(curElem,:);
	curP = x(curEdges);
	curCoords = c4n(n4e(curElem,:),:)';
	area = area4e(curElem);
	
    grad4e(curElem,:) = curP;
    
	signum = ones(3,1);
	I = find(e4ed(curEdges,2) == curElem);
	signum(I) = -1;
		
	psi1P1 = signum(1)*length4ed(curEdges(1))/(2*area)*(curCoords(:,1)-curCoords(:,3));	
	psi3P1 = signum(3)*length4ed(curEdges(3))/(2*area)*(curCoords(:,1)-curCoords(:,2));
	
	psi1P2 = signum(1)*length4ed(curEdges(1))/(2*area)*(curCoords(:,2)-curCoords(:,3));	
	psi2P2 = signum(2)*length4ed(curEdges(2))/(2*area)*(curCoords(:,2)-curCoords(:,1));
	
	psi2P3 = signum(2)*length4ed(curEdges(2))/(2*area)*(curCoords(:,3)-curCoords(:,1));	
	psi3P3 = signum(3)*length4ed(curEdges(3))/(2*area)*(curCoords(:,3)-curCoords(:,2));
	

	grad(1,:,curElem) = (curP(1)*psi1P1 + curP(3)*psi3P1);
	grad(2,:,curElem) = (curP(1)*psi1P2 + curP(2)*psi2P2);
	grad(3,:,curElem) = (curP(2)*psi2P3 + curP(3)*psi3P3);
	
end

%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).u = u;
p.level(end).u4e = u;
p.level(end).grad = grad;
p.level(end).grad4e = grad4e;
p.level(end).gradP = gradP;