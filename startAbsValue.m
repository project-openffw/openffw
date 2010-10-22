function p = startAbsValue
% start script for the nonlinear design problem
%
% author: Lena Noack 

%% discretization methods
pdeSolver = 'AVP1';

%% set up method specific configuration
% set up solver for the nonlinear system

% AVdirectFsolve solves the conforming method AVP1

solver = 'AVdirectFsolve';     %solves system of nonlinear equations F(x)=0
% solver = 'AVdirectFminunc';    %min of unconstrained multivariable function
% solver = 'AVdirectLineSearch'; %Newton-Raphson Scheme

% set up error estimator

% AVP1
%estimator = 'estimate_AvgS';    % ||sigma-A(sigma)||^1/2_L2, sigma=DW(grad Uh)
%estimator = 'estimate_Avg';    % ||ph-A(ph)||^1/2_L2, ph=grad Uh
estimator = 'estimate_Jump';   % |E|^1/2 [ph]_E^1/2
%estimator = 'estimate';        % |E| [grad Uh]_E + osc
%estimator = 'estimate2';       % ( h_T^2 ||f||^2_L^2 + 
                                %1/2 \sum_E(T) |E| ||[ph]_E\nu||^2_L^2 )^1/2

%% choose the problem
% numerical test example
problem = 'AbsValue_Square_exact';
% problem = 'AbsValue_Lshape_exact';
% problem = 'AbsValue_TwoSquare';

% set up parameters
maxNrDoF = 1000;
minDoF = 10;

% set up marking strategy
% mark = 'uniform'
mark = 'bulk'
% mark = 'max'

% Degree of Absolute Value: W(x)=|x|^p + regParam * |x|^2
pv = 3;
regParam = 0.1; %0.1;
fprintf('Degree of Absolute Value: W(x)=|x|^%.2g + %.5g*|x|^2 \n',pv,regParam);

% parameters for the nonlinear solver 'fsolve' of MATLAB
options = optimset('Display','off','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',50,'TolFun',1e-10,'TolX',1e-10);
%options = optimset('Display','iter','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',10,'TolFun',1e-10,'TolX',1e-10);

%% COMPUTE DISCRETE SOLUTION
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'AbsValue',solver,'redGreenBlue',estimator);
p.params.options = options;
p.params.RHS = 'RHS2';
p.problem.pv=pv;
p.problem.regParam=regParam;


dummy = (0.01 : 0.125 : 1.01).^3;
p.params.epsSequence = dummy(end:-1:1);
p.params.epsPower = 1;
p.params.epsFactor = 1;
p.params.rhsIntegtrateExactDegree = 19;
p.params.nonLinearExactIntegrateDegree = 5;
p = p.statics.run(p);

%% visualizing/post processing of the solution
 hold all

 figure(1);
 p = show('drawGrid',p);

 p.params.output.fontSize = 10;
 p.params.output.holdIt = false;
 p.params.output.saveFigures = true;
 p.params.output.minDoF = minDoF;
% p = drawConvergence(p)
% p = drawAllConv(p)
 p = drawFigures(p);
 p = evalFigures(p);

 
 
 
 

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

%figure(8);
%p = show('drawU',p);  %Discrete P1-solution bei ODRTGOAL
figure(9);
p = show('drawUexact',p);  %Discrete P1-solution bei ODRTGOAL
figure(10);
p = show('AVP1drawUminusUh',p);  %Discrete P1-solution bei ODRTGOAL

%figure(11);
%p = show('drawGradU',p); %Pfeile -> nur mit wenigen Freiheitsgraden %existiert noch nicht

% cd('results');
% mkdir('Test12');
% cd('Test12');
 
% save 'TwoSquare_ODP2_estEdge_h-dep_estimate_2500_bulk' p;
% cd('..')
% cd('..')

function p = evalFigures(p)

hold all
figure(13);
p.params.output.name = 'estError, uniform with \eta_{AvgSigma}';
p = show('drawError_estimatedError',p);
hold all
p.params.output.name = '||\sigma_h-\sigma||_{L^2(\Omega)}, uniform with \eta_{AvgSigma}';
p = show('drawError_L2errorPhminusP0',p);

grid off


%drawTriangleNew('Square',0.5,1)
