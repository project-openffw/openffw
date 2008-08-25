function teachRefinement
% visualisation of the mesh refinement

% Copyright 2007 Joscha Gedicke, Andreas Byfut
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


problem = 'Elliptic_Lshape_exact';
pdeSolver = 'P1';
maxNrDoF = 1000;

mark = 'maximum';

maxLevel = 6;

p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');
p.params.maxLevel = maxLevel;
p = computeSolution(p);

nrLvl = p.level(end).level;
p.params.output.holdIt = true;

for curLvl = 1 : nrLvl-1
   clf
   p.params.output.lineWidth = 2;
   drawGrid(p,curLvl);
   p.params.output.lineWidth = 2;
   drawRefEdges(p,curLvl); pause
   p.params.output.lineWidth = 4;
   drawMarkedBeforeClosureGrid(p,curLvl); pause
   p.params.output.lineWidth = 4;
   drawMarkedGrid(p,curLvl);pause
   p.params.output.lineWidth = 2;
   drawRGB(p,curLvl);pause
   p.params.output.lineWidth = 2;
   drawGrid(p,curLvl+1); pause
end
