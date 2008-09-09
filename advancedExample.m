function waterFallExample
% This script illustrates the performance of adaptive algorithms for an
% example that exhibits a strong gradient on a curve in the interior of
% the domain. Note that the strong gradient, where the magnitude can be
% controlled by the parameter k, does not result in a singular solution,
% i.e. uniform refinement results in an optimal convergence rate. Still,
% the adaptive algortihm gives better results with fewer degrees of freedom.

% The function we consider has the form:
% u = x*y*(1-x)*(1-y)*atan(k*(sqrt((x-5/4)^2 + (y+1/4)^2)-1));

k = 5;
computeAndPlot('uniform',k);
computeAndPlot('bulk',k);
% 
k = 50;
computeAndPlot('uniform',k);
computeAndPlot('bulk',k);


%%
function p = computeAndPlot(mark,k)

problem = 'Elliptic_Waterfall_exact';
pdeSolver = 'P1';         
maxNrDoF = 900;
integrationAccuracy = 19;

%% computes the solution
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elliptic');
p.params.rhsIntegtrateExactDegree = integrationAccuracy;
p.PDE.k = k;
p = p.statics.run(p);

%% plots the results in 3 figures: 
p.params.output.name = [mark,', ','k = ',num2str(k)];
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
