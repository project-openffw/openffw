function p = startTestreihe

% pdeSolver = 'ODP1' 'ODP2' 'ODRT' 'ODRTGOAL' 

% solver = 'ODlimitEps' 'ODoptimalLambda' 'ODdirectFsolve'
% 'ODeps_h_dependence'

% ODP1  % estimator = 'estimate_Avg' 'estimate_Jump';
% ODP2  % estimator = 'estimate' 'estimate2' 'estimateA'
% ODRT  % estimator = 'estimate_Proj' 'estimate_Jump';
% ODRTGOAL % estimator = 'estimate';

% problem = 'OptimalDesign_Square' 'OptimalDesign_Lshape' ;
% problem = 'OptimalDesign_Octagon' 'OptimalDesign_SquareSlit';
% problem = 'OptimalDesign_Square_exact' 'OptimalDesign_TwoSquare';

%Square % lambda = 0.0084;
%L-shape % lambda = 0.0143;
%Octagon % lambda = 0.0284;
%Square-Slit %lambda = 0.0168;

% Lagrange mupltiplier for numerical test example
% Square_Exact
% lambda = 0.003101917488281;

maxNrDoF = 3200;

% mark = 'uniform' 
mark = 'bulk'

mu1 = 1;
mu2 = 100;


% parameters for the nonlinear solver 'fsolve' of MATLAB
options = optimset('Display','iter','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',10,'TolFun',1e-10,'TolX',1e-10);
 
 
 
 
%% COMPUTE DISCRETE SOLUTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdeSolver = 'ODP2';
solver = 'ODeps_h_dependence';
estimator = 'estimate';
problem = 'OptimalDesign_TwoSquare';
lambda = 0.003101917488281;

t1 = (2*lambda*mu1/mu2).^(1/2);
t2 = mu2/mu1*t1;

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
p = p.statics.run(p);


%% visualizing/post processing of the solution
 hold all

 p.params.output.saveFigures = true;


 figure(1);
 p.params.output.name = 'estimated';
 p = show('drawError_estimatedError',p); 
 figure(2);
 p.params.output.name = 'L2';
 p = show('drawError_L2errorDisplacement',p);
 figure(3);
 p.params.output.name = 'H1semi';
 p = show('drawError_L2errorGradU',p);
 figure(4);
 p.params.output.name = 'NonLinear';
 p = show('drawError_L2errorPhminusP0',p);
 figure(5);
 p.params.output.name = 'Energy';
 p = show('drawError_errorEnergy',p);


 figure(6);
 p = show('drawGrid',p);
 figure(7);
 p = show('drawU',p);  %Discrete P1/P2-solution bei ODRTGOAL
 figure(8);
 p = show('ODP2drawUminusUh',p);  %Discrete P2-solution bei ODRTGOAL

% mkdir('results');
 cd('results');
 mkdir('Test10');
 cd('Test10');
 
 save 'TwoSquare_ODP2_h-dep_estimate_3000_bulk' p;
 cd('..')
 cd('..')


 
 
 
 
 
 
 mark = 'uniform'
 
 % parameters for the nonlinear solver 'fsolve' of MATLAB
options = optimset('Display','iter','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',10,'TolFun',1e-10,'TolX',1e-10);

%% COMPUTE DISCRETE SOLUTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdeSolver = 'ODP2';
solver = 'ODeps_h_dependence';
estimator = 'estimate';
problem = 'OptimalDesign_TwoSquare';
lambda = 0.003101917488281;

t1 = (2*lambda*mu1/mu2).^(1/2);
t2 = mu2/mu1*t1;

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
p = p.statics.run(p);


%% visualizing/post processing of the solution
 hold all

 p.params.output.saveFigures = true;


 figure(1);
 p.params.output.name = 'estimated';
 p = show('drawError_estimatedError',p); 
 figure(2);
 p.params.output.name = 'L2';
 p = show('drawError_L2errorDisplacement',p);
 figure(3);
 p.params.output.name = 'H1semi';
 p = show('drawError_L2errorGradU',p);
 figure(4);
 p.params.output.name = 'NonLinear';
 p = show('drawError_L2errorPhminusP0',p);
 figure(5);
 p.params.output.name = 'Energy';
 p = show('drawError_errorEnergy',p);


 figure(9);
 p = show('drawGrid',p);
 figure(10);
 p = show('drawU',p);  %Discrete P1/P2-solution bei ODRTGOAL
 figure(11);
 p = show('ODP2drawUminusUh',p);  %Discrete P2-solution bei ODRTGOAL

% mkdir('results');
 cd('results');
% mkdir('Test10');
 cd('Test10');
 
 save 'TwoSquare_ODP2_h-dep_estimate_3000_uniform' p;
 cd('..')
 cd('..')

