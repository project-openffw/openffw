function p = startElasticity
% start script for all elasticity problems

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


%% Problem Type
% problem = 'Elasticity_Square_exact';
% problem = 'Elasticity_Square_Neumann_exact';
problem = 'Elasticity_Lshape_exact';
% problem = 'Elasticity_Cooks';

%% PDE Solver
pdeSolver = 'P1P1';        %'P1-P1-Elasticity';
% pdeSolver = 'P1CR';       %'P1-CR-Elasticity';
% pdeSolver = 'AW';         %'Arnold-Winther-mixed-Elasticity';

%% Maximum Number Degrees of Freedom
maxNrDoF = 1000;

%% Mark Criterion
mark = 'bulk';
% mark = 'maximum';
% mark = 'uniform';
% mark = 'graded';

%% Compute Discrete Solution

p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elasticity');
% Overwrites the standard parameters (nu,E) defined in the problem definition
p = elasticity_setLame(p,0.3,1000);
p = computeSolution(p);

%% Show Discrete Solution
figure(1)
p.params.output.holdIt = false;
p.params.output.lineWidth = 1;
p.params.output.factor = 5000;

set(gcf,'Name','Displacement');
p = show('drawU',p);

%% Show Grid
figure(2)
p.params.output.holdIt = false; 

set(gcf,'Name','Grid');
p = show('drawGrid',p);

%% Show Convergence History
figure(3)
set(gcf,'Name','Convergence History');
p.params.output.drawConvergenceRate = false;
% p.params.integrationDegrees.exactError = 9;
p.params.output.holdIt = true; 

p.params.output.name = [mark,', ',pdeSolver,', ','estimator'];
p = show('drawError_estimatedError',p);
if ~isempty(findstr(problem,'exact'))
    p.params.output.name = [mark,', ',pdeSolver,', ','H^1'];
    p = show('drawError_H1semiError',p);
    p.params.output.name = [mark,', ',pdeSolver,', ','L^2'];
    p = show('drawError_L2error',p);
end