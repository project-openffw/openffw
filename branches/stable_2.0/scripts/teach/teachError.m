function teachError
% draws the errror on the grid

% Copyright 2007 Andreas Byfut
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
mark = 'bulk';
maxNrDoF = 1000;

maxLevel = 10; 

p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');
p.params.maxLevel = maxLevel;
p = computeSolution(p);

p.params.errorIntegtrateExactDegree = 9;
p.params.IntegrationMode = 'vectorised';

nrLvl = p.level(end).level;
p.params.output.holdIt = false;

for curLvl = 1 : nrLvl
   drawErrorOnGrid(p,curLvl,'H1semiError');
   view(-30,20)
   pause   
end
