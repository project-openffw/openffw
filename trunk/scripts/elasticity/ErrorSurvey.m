function ErrorSurvey

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
% problem = 'Elasticity_Lshape_exact';
problem = 'Elasticity_Square_exact';

%% Degrees of Freedom
maxNrDoF = 501;

%% mark
% mark = 'bulk';	
mark = 'uniform';

%% Parameters
nu = 0.3;

exactH1semiErrorTolerance = 1;
exactL2errorTolerance = 1e-8;

%% Error Type
errorType = 'drawError_L2error';
% errorType = 'drawError_H1semiError';
% errorType = 'drawError_estimatedError';


%% Solve

ErrorComparison(problem,mark,maxNrDoF,nu,errorType,...
    exactL2errorTolerance,exactH1semiErrorTolerance);


%%
function ErrorComparison(problem,mark,maxNrDoF,nu,errType,...
        exactL2errorTolerance,exactH1semiErrorTolerance)



%% COMPUTE DISCRETE SOLUTION
pdeSolver = 'AW';
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elasticity');
p.params.exactL2errorTolerance = exactL2errorTolerance;
p.params.exactH1semiErrorTolerance = exactH1semiErrorTolerance;
p = elasticity_setLame(p,nu);
p = computeSolution(p);

p.params.output.holdIt = true;
p.params.output.plotGrid = false;
p.params.output.name = 'AW';
p.params.output.myColor = 'r';
show(errType,p);


%% COMPUTE DISCRETE SOLUTION
pdeSolver = 'P1P1';
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elasticity');
p.params.exactL2errorTolerance = exactL2errorTolerance;
p.params.exactH1semiErrorTolerance = exactH1semiErrorTolerance;
p = elasticity_setLame(p,nu);
p = computeSolution(p);

p.params.output.holdIt = true;
p.params.output.plotGrid = false;
p.params.output.name = 'P1P1';
p.params.output.myColor = 'g';
show(errType,p);

%% COMPUTE DISCRETE SOLUTION
pdeSolver = 'P1CR';
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elasticity');
p.params.exactL2errorTolerance = exactL2errorTolerance;
p.params.exactH1semiErrorTolerance = exactH1semiErrorTolerance;
p = elasticity_setLame(p,nu);
p = computeSolution(p);

p.params.output.holdIt = true;
p.params.output.plotGrid = false;
p.params.output.name = 'P1CR';
p.params.output.myColor = 'b';
show(errType,p);



