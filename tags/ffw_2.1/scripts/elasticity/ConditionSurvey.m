function ConditionSurvey

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
problem = 'Elasticity_Lshape_exact';
% problem = 'Elasticity_Square_exact';

%% Degrees of Freedom
maxNrDoF = 1400;

%% Mark
% mark = 'bulk';	
mark = 'uniform';

%% nu
nu1 = 0.3;
nu2 = 0.49;
nu3 = 0.4999;

%% Solve
clf
CondComparison('AW',problem,mark,maxNrDoF,nu1,nu2,nu3,'x');
CondComparison('P1CR',problem,mark,maxNrDoF,nu1,nu2,nu3,'o');
CondComparison('P1P1',problem,mark,maxNrDoF,nu1,nu2,nu3,'*');



%%
function CondComparison(pdeSolver,problem,mark,maxNrDoF,nu1,nu2,nu3,marker)

%% COMPUTE DISCRETE SOLUTION
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elasticity');
p.params.showIteratively = false;
p = elasticity_setLame(p,nu1);
p = computeSolution(p);

p.params.output.holdIt = true;
p.params.output.plotGrid = false;
p.params.output.name = ['nu = ',num2str(nu1),', ',pdeSolver];
p.params.output.myColor = 'r';
p.params.output.marker = marker;
show('drawCondNr',p);


%% COMPUTE DISCRETE SOLUTION
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elasticity');
p.params.showIteratively = false;
p = elasticity_setLame(p,nu2);
p = computeSolution(p);

p.params.output.holdIt = true;
p.params.output.plotGrid = false;
p.params.output.name = ['nu = ',num2str(nu2),', ',pdeSolver];
p.params.output.myColor = 'g';
p.params.output.marker = marker;
show('drawCondNr',p);

%% COMPUTE DISCRETE SOLUTION
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'elasticity');
p.params.showIteratively = false;
p = elasticity_setLame(p,nu3);
p = computeSolution(p);

p.params.output.holdIt = true;
p.params.output.plotGrid = false;
p.params.output.name = ['nu = ',num2str(nu3),', ',pdeSolver];
p.params.output.myColor = 'b';
p.params.output.marker = marker;
show('drawCondNr',p);


