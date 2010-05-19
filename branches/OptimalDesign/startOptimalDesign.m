function p = startOptimalDesign
% Copyright 2007 David Guenther
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
%
%
%
% start script for the nonlinear design problem to find the optimal
% material distribution
%
% Three discritization methods are available: the conforming P1- and P2-method 
% and the mixed RT0xP0 method. 
%
% Besides the discritization methods and a direct
% solve of the problem you can also choose to solve the
% problem by taking a cost functional J which has to be minimize under the
% restriction of the weak formula of the problem.
%
% One can choose between a uniform or an adaptive refinement. Also there
% are several algorithms/solver available to compute the solution.
%
% author: David Guenther 

%% discretization methods
% pdeSolver = 'ODP1';
pdeSolver = 'ODP2';
% pdeSolver = 'ODRT';

%% goal oriented discretization methods
% minimize the cost functional:    J(sigma) = \int W*(|\sigma|) dx
% pdeSolver = 'ODRTGOAL';

%% set up method specific configuration
% set up solver for the nonlinear system

% the first two solve the regularized system for the mixed formulation, 
% in ODeps_h_dependence the regularization parameter is mesh-dependent, in 
% ODlimitEps the limit of the regularization is approximated by a sequence 
% of regularization parameters; ODoptimalLambda determine the 
% problem-specific optimal Lagrange-Multiplier for the design constraint 
% (use uniform refinement only); 
% ODdirectFsolve solves the conforming methods like ODP1 and ODP2

% solver = 'ODeps_h_dependence';
% solver = 'ODlimitEps';
% solver = 'ODoptimalLambda';
solver = 'ODdirectFsolve';

% set up error estimator

% For the conforming P1 method you can choose between an average estimator
% and the estimator defined by the the jump of stress; for the P2 method
% just the jump estimator is available; in the mixed formulation one can
% choose between the projection estimator (computation of the
% rotationfield) and a jump estimator

% ODP1
% estimator = 'estimate_Avg';
% estimator = 'estimate_Jump';
% ODP2 / ODRTGOAL
estimator = 'estimate';
% ODRT
% estimator = 'estimate_Proj';
% estimator = 'estimate_Jump';

%% choose the problem
% model problems
% problem = 'OptimalDesign_Square';
problem = 'OptimalDesign_Lshape';
% problem = 'OptimalDesign_Octagon';
% problem = 'OptimalDesign_SquareSlit';

% numerical test example
% problem = 'OptimalDesign_Square_exact';

% set up parameters
maxNrDoF = 1000;

% set up marking strategy
% mark = 'uniform'
mark = 'bulk'
% mark = 'max'

% material parameters mu_1 and mu_2
mu1 = 1;
mu2 = 2;

% Lagrange mupltiplier for the model problems
%Square
% lambda = 0.0084;
%L-shape
lambda = 0.0143;
%Octagon
% lambda = 0.0284;
%Square-Slit
% lambda = 0.0168;

% Lagrange mupltiplier for numerical test example
% Square_Exact
% lambda = 0.003101917488281;

t1 = (2*lambda*mu1/mu2).^(1/2);
t2 = mu2/mu1*t1;

% parameters for the nonlinear solver 'fsolve' of MATLAB
options = optimset('Display','iter','Jacobian','on','NonlEqnAlgorithm','dogleg','MaxIter',10,'TolFun',1e-10,'TolX',1e-10);

%% COMPUTE DISCRETE SOLUTION
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'optimalDesign',solver,'redGreenBlue',estimator);
p.params.options = options;
p.params.RHS = 'RHS2';
p.problem.mu1 = mu1;
p.problem.mu2 = mu2;

p.problem.t1 = t1;
p.problem.t2 = t2;
p.problem.lambda = lambda;

dummy = (0.01 : 0.125 : 1.01).^3;
p.params.epsSequence = dummy(end:-1:1);
p.params.epsPower = 1;
p.params.epsFactor = 1;
p.params.rhsIntegtrateExactDegree = 19;
p.params.nonLinearExactIntegrateDegree = 5;
p = ODcomputeSolution(p);

%% visualizing/post processing of the solution
% hold all

% p = show('drawGrid',p);
p = ODshowVolumeFraction(p,p.level(end).level);

p.params.output.fontSize = 12;
p.params.output.holdIt = true;
% p = drawFigures(p);

%% draw convergence history
function p = drawFigures(p)

% set(gcf,'Name','Square exact 2');
p.params.output.name = 'estimated';
p = show('drawError_estimatedError',p);
p.params.output.name = 'L2';
p = show('drawError_L2errorDisplacement',p);
p.params.output.name = 'H1semi';
p = show('drawError_L2errorGradU',p);
p.params.output.name = 'NonLinear';
p = show('drawError_L2errorPhminusP0',p);
p.params.output.name = 'Energy';
p = show('drawError_errorEnergy',p);

