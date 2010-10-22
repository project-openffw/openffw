function p = startTwoWell
% start script for the nonlinear design problem
%
% author: Lena Noack 

%% discretization methods
pdeSolver = 'TWP1';
%pdeSolver = 'TWP2';
%pdeSolver = 'TWDWRP1';
%pdeSolver = 'TWDWR';
%pdeSolver = 'TWRT';

%% set up method specific configuration
% set up solver for the nonlinear system

% the first two solve the regularized system for the mixed formulation, 
% in TWeps_h_dependence the regularization parameter is mesh-dependent, in 
% TWlimitEps the limit of the regularization is approximated by a sequence 
% of regularization parameters; TWoptimalLambda determine the 
% problem-specific optimal Lagrange-Multiplier for the design constraint; 
% TWdirectFsolve solves the conforming method TWP1

% solver = 'TWeps_h_dependence';
% solver = 'TWlimitEps';
% solver = 'TWoptimalLambda';
solver = 'TWdirectFsolve';
% solver = 'TWdirectLineSearch'; %Newton-Raphson Scheme

% set up error estimator

% P1
% estimator = 'estimate';
% estimator = 'estimate2';
% estimator = 'estimate_Avg';
% estimator = 'estimate_Jump';
 estimator = 'estimate_Jump_SigmaP'; %ok 
% estimator = 'estimate_Jump_SigmaPL2'; %ok 
% estimator = 'estimate_Jump_Avg'; %ok bei OD, hier nicht
% estimator = 'estimate_JumpMark';
% estimator = 'estimateEdOsc';

% P2
% estimator = 'estimate';
% estimator = 'estimate_Jump_SigmaP';
% estimator = 'estimate_Jump_SigmaPL2';

% DWR-P1
% estimator = 'estimateDWR';

% P0-RT
 estimator = 'estimate';
% estimator = 'estimate_Jump_SigmaP';

%% choose the problem
% numerical test example
problem = 'TwoWell_Square';
% problem = 'TwoWell_TwoSquare';
% problem = 'TwoWell_Lshape_exact';
% problem = 'TwoWell_Lshape';
% problem = 'TwoWell_Octagon';
% problem = 'TwoWell_SquareSlit';

% set up parameters
maxNrDoF = 100;
minDoF = 10;

% set up marking strategy
% mark = 'uniform'
 mark = 'bulk'
% mark = 'max'

% Boundaries F_1 and F_2;
% F_1 < F_2 < infty
F1 = -(1/sqrt(13))*[3;2];
F2 = (1/sqrt(13))*[3;2];

% parameters for the nonlinear solver 'fsolve' of MATLAB
options = optimset('Display','off','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',10,'TolFun',1e-10,'TolX',1e-10);
%options = optimset('Display','off','Jacobian','on','NonlEqnAlgorithm','lm','MaxIter',30,'TolFun',1e-10,'TolX',1e-10);
%options = optimset('Display','iter','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',10,'TolFun',1e-10,'TolX',1e-10);

%% COMPUTE DISCRETE SOLUTION
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'TwoWell',solver,'redGreenBlue',estimator);
p.params.options = options;
p.params.RHS = 'RHS2';
p.params.CONV = 'c'; %convex energy density
%p.params.CONV = 'nc'; %non-convex energy density
p.problem.F1 = F1;
p.problem.F2 = F2;

disp([num2str(maxNrDoF) ' Freiheitsgrade, W=' p.params.CONV]);

%CONV = p.params.CONV;
%if strcmp(CONV,'c')
%else
%end
        
dummy = (0.01 : 0.125 : 1.01).^3;
p.params.epsSequence = dummy(end:-1:1);
p.params.epsPower = 1;
p.params.epsFactor = 1;
p.params.rhsIntegtrateExactDegree = 19;
p.params.nonLinearExactIntegrateDegree = 5;
p = p.statics.run(p);

%% visualizing/post processing of the solution
 hold all

 p = show('drawGrid',p);

 p.params.output.fontSize = 10;
 p.params.output.holdIt = false;
 p.params.output.saveFigures = true;
 p.params.output.minDoF = minDoF;
% p = drawConvergence(p)
 p = drawAllConv(p)
 p = drawFigures(p);
% p = evalFigures(p);
 
 
 
 

%% draw convergence history
function p = drawConvergence(p)

 set(gcf,'Name','Square');
figure(2);
p.params.output.name = 'estError mit \eta_{Jump}';
p = show('drawError_estimatedError',p);
figure(3);
p.params.output.name = '||u-u_h||_{L^2(\Omega)}';
p = show('drawError_L2errorDisplacement',p);
figure(4);
p.params.output.name = '||\nabla u-\nabla u_h||_{L^2(\Omega)}';
p = show('drawError_L2errorGradU',p);
figure(5);
p.params.output.name = '||\sigma_h-\sigma||_{L^2(\Omega)}';
p = show('drawError_L2errorPhminusP0',p);
figure(6);
p.params.output.name = '|E(u_h)-E(u)|';
p = show('drawError_errorEnergy',p);

grid off

function p = drawAllConv(p)

figure(12);
p.params.output.name = 'estError mit \eta^{(1)}';
%p.params.output.name = 'estError mit \eta_{Jump}';
p = show('drawError_estimatedError',p);
hold all
p.params.output.name = '||u-u_h||_{L^2(\Omega)}';
p = show('drawError_L2errorDisplacement',p);
hold all
p.params.output.name = '||\nabla u-\nabla u_h||_{L^2(\Omega)}';
p = show('drawError_L2errorGradU',p);
hold all
p.params.output.name = '||\sigma_h-\sigma||_{L^2(\Omega)}';
p = show('drawError_L2errorPhminusP0',p);
hold all
p.params.output.name = '|E(u_h)-E(u)|';
p = show('drawError_errorEnergy',p);

grid off

function p = drawFigures(p)

figure(8);
p = show('drawU',p);  %Discrete solution bei ODRTGOAL
figure(9);
p = show('drawUexact',p);  %Discrete solution bei ODRTGOAL
figure(10);
p = show('drawUminusUh',p);  %Discrete solution bei ODRTGOAL

%figure(11);
%p = show('drawGradU',p); %Pfeile -> nur mit wenigen Freiheitsgraden %existiert noch nicht

% cd('results');
% mkdir('Test12');
% cd('Test12');
 
% save 'TwoSquare_ODP2_estEdge_h-dep_estimate_2500_bulk' p;
% cd('..')
% cd('..')

function p = evalFigures(p)

figure(13);
p.params.output.name = 'estError mit \eta_{Jump}';
p = show('drawError_estimatedError',p);
hold all
figure(14);
p.params.output.name = '||\sigma_h-\sigma||_{L^2(\Omega)} mit \eta_{Jump}';
p = show('drawError_L2errorPhminusP0',p);

grid off

%drawTriangleNew('Square',0.5,1)

