function p = example5

% Copyright 2007 Joscha Gedicke
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


%% Problem Type
problem = 'Elliptic_Waterfall_exact';

%% PDE Solver
pdeSolver = 'P1';

%% Maximum Number Degrees of Freedom
maxNrDoF = 50;

%% Mark Criterion
mark = 'bulk';

%% Compute Discrete Solution
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');

p.PDE.k = 120;  %<-- change this parameter and observe the 
               %    change in the convergence history

p = computeSolution(p);

%% Show Discrete Solution
figure
p.params.output.holdIt = false; 
   % % Due to problems when plotting 3D graphics of large DoF
   % prntlvl=1;
   % for i=1:size(p.level,2)
   % 	   if p.level(i).nrElems<200
   % 		   prntlvl=i;
   % 	   else
   % 		   break;
   % 	   end
   % end
   % p = show('drawU',p,prntlvl);
   % Otherwise, to plot the solution of the last level
   p = show('drawU',p);

%% Show Grid
figure
p.params.output.holdIt = false; 
p = show('drawGrid',p);

%% Show Convergence History
figure

p.params.output.minDoF = 1;
p.params.output.drawConvergenceRate = false;
p.params.errorIntegtrateExactDegree = 9;
p.params.output.holdIt = true; 

p.params.output.name = [mark,', ',pdeSolver,', ','estimator'];
p = show('drawError_estimatedError',p);

p.params.IntegrationMode = 'vectorised';
p.params.output.name = [mark,', ',pdeSolver,', ','H^1'];
p = show('drawError_H1semiError',p);
p.params.output.name = [mark,', ',pdeSolver,', ','L^2'];
p = show('drawError_L2error',p);


