function advancedExample
% This script illustrates the performance of adaptive algorithms for an
% example that exhibits a strong gradient on a curve in the interior of
% the domain. Note that the strong gradient, where the magnitude can be
% controlled by the parameter k, does not result in a singular solution,
% i.e. uniform refinement results in an optimal convergence rate. Still,
% the adaptive algortihm gives better results with fewer degrees of freedom.
%
% The function we consider has the form:
% u = x*y*(1-x)*(1-y)*atan(k*(sqrt((x-5/4)^2 + (y+1/4)^2)-1));

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


K = 5;
computeAndPlot('uniform',K);
computeAndPlot('bulk',K);
% 
K = 50;
computeAndPlot('uniform',K);
computeAndPlot('bulk',K);


%%
function p = computeAndPlot(mark,K)

problem = 'Elliptic_Waterfall_exact';
pdeSolver = 'P1';         
maxNrDoF = 900;
integrationAccuracy = 9;

%% computes the solution
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');
p.params.rhsIntegtrateExactDegree = integrationAccuracy;
p.params.errorIntegtrateExactDegree = integrationAccuracy;
p.PDE.K = K;

p = computeSolution(p);

%% plots the results in 3 figures: 
p.params.output.name = [mark,', ','K = ',num2str(K)];
p.params.output.holdIt = true; 

figure(1)
set(gcf,'Name','Displacement');
p = show('drawU',p);

figure(2)
set(gcf,'Name','Displacement Error');
p = show('drawError_L2error',p);

figure(3)
set(gcf,'Name','Energy Error');
p = show('drawError_H1semiError',p);
