function p = startOptimalDesign
% start script for the nonlinear design problem to find the optimal
% material distribution
%
% Three discritization methods are available: the conforming P1- and P2-method 
% and the mixed RT0xP0 method. 
%
% Besides the discritization methods and a direct
% solve of the problem you can also choose to solve the
% problem by taking a cost functional J which has to be minimized under the
% restriction of the weak formula of the problem.
%
% One can choose between a uniform or an adaptive refinement. Also there
% are several algorithms/solver available to compute the solution.
%
% author: David Guenther 

%% discretization methods
pdeSolver = 'ODP1';
% pdeSolver = 'ODP2';
% pdeSolver = 'ODRT';

%% goal oriented discretization methods
% minimize the cost functional:    J(sigma)
% pdeSolver = 'ODRTGOAL';

%% set up method specific configuration
% set up solver for the nonlinear system

% the first two solve the regularized system for the mixed formulation, 
% in ODeps_h_dependence the regularization parameter is mesh-dependent, in 
% ODlimitEps the limit of the regularization is approximated by a sequence 
% of regularization parameters; ODoptimalLambda determine the 
% problem-specific optimal Lagrange-Multiplier for the design constraint; 
% ODdirectFsolve solves the conforming methods like ODP1 and ODP2

% solver = 'ODeps_h_dependence';
% solver = 'ODlimitEps';
% solver = 'ODoptimalLambda';
solver = 'ODoptimalLambda2';
% solver = 'ODoptimalLambdaSeq';
% solver = 'ODdirectFsolve';
% solver = 'ODdirectFminunc'; %funkt. nicht
% solver = 'ODdirectLineSearch';

% set up error estimator

% For the conforming P1 method you can choose between an average estimator
% and the estimator defined by the the jump of stress; for the P2 method
% just the jump estimator is available; in the mixed formulation one can
% choose between the projection estimator (computation of the
% rotationfield) and a jump estimator

% ODP1
% estimator = 'estimateInterp';
estimator = 'estimateH2';
% estimator = 'estimateCompareSols';
% estimator = 'estimateDeltaEnergy';

% estimator = 'estimate_PD';
% estimator = 'estimate_Avg';
% estimator = 'estimate_Jump';
% estimator = 'estimate_JumpSquare';
% ODP2
% estimator = 'estimate';
% estimator = 'estimate2';
% estimator = 'estimateA';
% estimator = 'estimateEdge';
% ODRT
% estimator = 'estimate';
% estimator = 'estimate_Proj';
% estimator = 'estimate_Jump';
% ODRTGOAL
% estimator = 'estimate';
% estimator = 'estimate2';
% estimator = 'estimateA';

% New Error Estimators
%estimator = 'estimate_Jump_SigmaP'; %ok mit 1 h-Potenz weniger
%estimator = 'estimate_Jump_Avg'; %ok
%estimator = 'estimate_RHS_Jump'; %Noch Ju_h
%estimator = 'estimate_PSigma_lRHS'; %Funktioniert nicht
%estimator = 'estimate_PSigma_nlRHS'; %Funktioniert nicht

%% choose the problem
% model problems
% problem = 'OptimalDesign_Square';
% problem = 'OptimalDesign_Lshape';
% problem = 'OptimalDesign_Octagon';
% problem = 'OptimalDesign_SquareSlit';

% numerical test example
% problem = 'OptimalDesign_TwoSquare';
% problem = 'OptimalDesign_Square_exact';
% problem = 'OptimalDesign_Square_exact2';
% problem = 'OptimalDesign_Lshape_exact';
problem = 'OptimalDesign_SquareSlit_exact';

% set up parameters
maxNrDoF = 500;
minDoF = 5;

% set up marking strategy
mark = 'uniform'
% mark = 'bulk'
% mark = 'max'

% material parameters mu_1 and mu_2;
% 0 < mu_1 < mu_2 < infty
mu1 = 1;
mu2 = 2;

% Lagrange mupltiplier for the model problems
%Square
% lambda = 0.00845150868686 %nach Bartels/Carstensen, 4 x red
% energy =-0.01515320514492;
% lambda = 0.00844992886957; %mu1=1, mu2=2, 4 x red
%          0.00844993115756 mit P=0.9
% energy =-0.01515305599894;
%L-shape
% lambda = 0.01448967844544; %nach Bartels/Carstensen
% energy =-0.09527515596197;
% lambda = 0.01450891266029; %mu1=1, mu2=2, 4 x red
% energy =-0.09527515452131;
%Octagon
% lambda = 0.02844454033907; %nach Bartels/Carstensen
% energy =-0.13636758964217;
% lambda = 0.02128622918119; %mu1=1, mu2=2, 4 x red
% energy =-0.13758519182292;
%Square-Slit
% lambda = 0.01636820082193; %nach Bartels/Carstensen
% energy =-0.14336007558320;
% lambda = 0.01622726723795; %mu1=1, mu2=2, 4 x red
% energy =-0.14259590630143;

% Lagrange mupltiplier for numerical test examples
% TwoSquare
% lambda = 0; %mu1=1, mu2=2, 4 x red
% energy = 0;
% Square_Exact
% lambda = 0.003101917488281; %????
% lambda = 0.00331580978555; %mu1=1, mu2=2, 4 x red
% energy =-0.00999205270568;
% Square_Exact2
% lambda = 0.00582753505851; %mu1=1, mu2=2, 4 x red
% energy =-0.03245743852395;
% Lshape_Exact
% lambda = 0.36462636980728; %mu1=1, mu2=2, 4 x red
% energy = 1.80479706076594;
%Square-Slit
lambda = 0.013561899770084; %mu1=1, mu2=2, 4 x red
% energy =-0;


t1 = (2*lambda*mu1/mu2).^(1/2);
t2 = mu2/mu1*t1;

% parameters for the nonlinear solver 'fsolve' of MATLAB
%options = optimset('Display','iter','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',10,'TolFun',1e-10,'TolX',1e-10);
options = optimset('Display','off','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',50,'TolFun',1e-10,'TolX',1e-10);

%% COMPUTE DISCRETE SOLUTION
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'optimalDesign',solver,'redGreenBlue',estimator);
p.params.options = options;
p.params.RHS = 'RHS2';
p.problem.mu1 = mu1;
p.problem.mu2 = mu2;

p.problem.t1 = t1;
p.problem.t2 = t2;
p.problem.lambda = lambda;

dummy = (0.0 : 0.002 : 0.01).^3;
%dummy = (0.01 : 0.125 : 1.01).^3;
p.params.epsSequence = dummy(end:-1:1);
p.params.epsPower = 1;
p.params.epsFactor = 1;
p.params.rhsIntegtrateExactDegree = 19;
p.params.nonLinearExactIntegrateDegree = 5;
p = p.statics.run(p);


%% visualizing/post processing of the solution
 hold all

 figure(3)
 p = show('drawGrid',p);

 p.params.output.fontSize = 10;
 p.params.output.holdIt = false;
 p.params.output.saveFigures = true;
 p.params.output.minDoF = minDoF;
 
 % p = drawConvergence(p);
% p = drawAllConv(p);
p = drawFigures(p);
% p = evalFigures(p);
p = drawNoExact(p,estimator);

% EnMid = getMidEnergy(p)

%% draw convergence history
function p = drawConvergence(p)

 set(gcf,'Name','Square');
figure(2);
p.params.output.name = 'estError';
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
p.params.output.name = 'estError';
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
p = show('drawU',p);  %Discrete P1-solution bei ODRTGOAL
%figure(9);
%p = show('drawUexact',p);  %Discrete P1-solution bei ODRTGOAL
%figure(10);
%p = show('drawUminusUh',p);  %Discrete P1-solution bei ODRTGOAL



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
p.params.output.name = 'estError';
p = show('drawError_estimatedError',p);
hold all
figure(14);
p.params.output.name = '||\sigma_h-\sigma||_{L^2(\Omega)}';
p = show('drawError_L2errorPhminusP0',p);

grid off

function p = drawNoExact(p,estimator)


figure(1);
hold all
p.params.output.name = [estimator, '\delta, Lshape'];
p = show('drawError_estimatedError',p);

%figure(3);
%p = show('drawU',p);  %Discrete P1-solution bei ODRTGOAL

grid off

%drawTriangleNew('Square',0.5,1)

function EnMid = getMidEnergy(p)
length4ed = p.level(end).enum.length4ed;
nrElems = p.level(end).nrElems;
n4e = p.level(end).geom.n4e;
lvl = length(p.level);
ed4e = p.level(end).enum.ed4e;
h_T = max(length4ed(ed4e)')';
area4T = 0.25*h_T.^2;%p.level(end).enum.area4T;

intEnergy = integrate(n4e,lvl,10,@getEnergy,p);
EnTotal = sum(intEnergy);
EnMid = (1/nrElems)*sum(intEnergy./area4T);

fprintf('\nEnTotal = %.14g \t EnMid = %.14g\n',EnTotal, EnMid)
% Lshape
%EnTotal = -0.095276255198805 	 
%EnMid = -0.031758751732935


%% supply the discrete energy
function val = getEnergy(x,y,curElem,lvl,p)

energy_h = p.statics.energy_h;
val = energy_h(x,y,curElem,lvl,p);

