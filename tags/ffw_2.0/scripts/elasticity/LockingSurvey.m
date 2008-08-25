function LockingSurvey

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


%% Problem
problem = 'Elasticity_Square_exact';
% problem = 'Elasticity_Lshape_exact';
% problem = 'Elasticity_Cooks';

%% Degrees of Freedom
maxNrDoF = 1000;

%% Mark
mark = 'uniform';
% mark = 'bulk';


%% Error Type
errType = 'drawError_H1semiError';
% errType = 'drawError_L2Error';
% errType = 'drawError_estimatedError';

%% nu
nu1 = 0.3;
nu2 = 0.4999;
nu3 = 0.49999999;

%% Tolerance
exactH1semiErrorToleranceNU1 = 1;
exactH1semiErrorToleranceNU2 = 10;
exactH1semiErrorToleranceNU3 = 100;


%% Solve
pdeSolver = 'AW';
% pdeSolver = 'P1CR';
% pdeSolver = 'P1P1';

LockingComparison(pdeSolver,problem,mark,maxNrDoF,nu1,nu2,nu3,errType,...
     exactH1semiErrorToleranceNU1,exactH1semiErrorToleranceNU2,exactH1semiErrorToleranceNU3);
 

%%
function LockingComparison(pdeSolver,problem,mark,maxNrDoF,nu1,nu2,nu3,errType,...
             exactH1semiErrorToleranceNU1,exactH1semiErrorToleranceNU2,exactH1semiErrorToleranceNU3)



%% COMPUTE DISCRETE SOLUTION
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elasticity');
p.params.exactH1semiErrorTolerance = exactH1semiErrorToleranceNU1;
p.params.exactL2errorTolerance = exactH1semiErrorToleranceNU1;
p.params.showIteratively = false;
p = elasticity_setLame(p,nu1);
p = computeSolution(p);

p.params.output.holdIt = true;
p.params.output.plotGrid = false;
p.params.output.name = ['nu = ',num2str(nu1),', ',pdeSolver];
p.params.output.myColor = 'r';
p = show(errType,p);


%% COMPUTE DISCRETE SOLUTION
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elasticity');
p.params.exactH1semiErrorTolerance = exactH1semiErrorToleranceNU2;
p.params.exactL2errorTolerance = exactH1semiErrorToleranceNU2;
p.params.showIteratively = false;
p = elasticity_setLame(p,nu2);
p = computeSolution(p);

p.params.output.holdIt = true;
p.params.output.plotGrid = false;
p.params.output.name = ['nu = ',num2str(nu2),', ',pdeSolver];
p.params.output.myColor = 'g';
p = show(errType,p);

%% COMPUTE DISCRETE SOLUTION
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elasticity');
p.params.exactH1semiErrorTolerance = exactH1semiErrorToleranceNU3;
p.params.exactL2errorTolerance = exactH1semiErrorToleranceNU3;
p.params.showIteratively = false;
p = elasticity_setLame(p,nu3);
p = computeSolution(p);

p.params.output.holdIt = true;
p.params.output.plotGrid = false;
p.params.output.name = ['nu = ',num2str(nu3),', ',pdeSolver];
p.params.output.myColor = 'b';
p = show(errType,p);

