function p = startDWR(estimate,nr,DrawSigma,GOAL,lambda)
% Copyright 2008 Joscha Gedicke, Lena Noack

getLambda = 0;
BoolL43 = 0;   % 0 -> L^2-Norm, 1 -> L^4/3-Norm

if(nargin < 5)
    getLambda = 1;
end
if(nargin < 4)
    GOAL = 'GOAL_IM';
end
if(nargin < 3)
    DrawSigma = 0;
end
if(nargin < 2)
    nr = 3;
end
if(nargin < 1)
    estimate = 'estimateDeltaEnergy';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimal Design

%% Problem Type
%problem = 'OptimalDesign_SquareSlit';
%problem = 'OptimalDesign_SquareSlit_exact';
problem = 'OptimalDesign_Lshape';
%problem = 'OptimalDesign_Lshape_exact';

probClass = 'optimalDesign';

%% PDE Solver
pdeSolver = 'ODRT'; 
   DWRbool=1;
   %DWRbool=0;
%pdeSolver = 'ODP1'; 
   %DWRbool=1;
   %DWRbool=0;

% definition of constants for problem
mu1 = 1;
mu2 = 2;

if getLambda
  if strcmp(problem,'OptimalDesign_SquareSlit')
  %    lambda = 0.01622726723795;   % SquareSlit, P1
      lambda = 0.019549371741127;  % SquareSlit, PO-RT0 2x red
  elseif strcmp(problem,'OptimalDesign_SquareSlit_exact')
  %    lambda = 0.013561899770084;  % SquareSlit exact
      lambda = 0.025603456541462;  % SquareSlit exact, PO-RT0 2x red
  elseif strcmp(problem,'OptimalDesign_Lshape')
  %    lambda = 0.01622726723795;   % SquareSlit
  %    lambda = 0.012647410259876;  % Lshape, P1 1xred
  %    lambda = 0.015065000646679;  % Lshape, P1 2xred
  %    lambda = 0.01438421216343;   % Lshape, P1 3xred
  %    lambda = 0.014508912660285;  % Lshape, P1 4xred
  %    lambda = 0.0155;             % Lshape, P1 best results
      lambda = 0.013084509211941;  % Lshape, PO-RT0 2xred
  else 
  %    lambda = 0.36447407757502;  % Lshape exact, P1 1xred
  %    lambda = 0.36564158241509;  % Lshape exact, P1 2xred
  %    lambda = 0.36261385275783;  % Lshape exact, P1 3xred
  %    lambda = 0.3646263768783;   % Lshape exact, P1 4xred
      lambda = 0.36999663504332;  % Lshape exact, PO-RT0 2xred
  end
end

%% Method to find minimum of E(v)
solver = 'ODdirectFsolve';
%% or Method to find optimal Lambda
%solver = 'ODoptimalLambda2';

t1 = (2*lambda*mu1/mu2).^(1/2);
t2 = mu2/mu1*t1;

F1=0;F2=0;CONV = ''; pv=0; regParam=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2-Well Benchmark

%% Problem Type
%problem = 'TwoWell_Square'; %TWCF
%problem = 'TwoWell_SquareExact'; %TWCO

%probClass = 'TwoWell'; BoolL43 = 1;

%% PDE Solver
%pdeSolver = 'TWP1';    DWRbool=0;
%pdeSolver = 'TWP1DWR'; DWRbool=1;
%pdeSolver = 'TWP2';    DWRbool=0;

% definition of constants for problem
% Boundaries F_1 and F_2;
% F_1 < F_2 < infty
%F1 = -(1/sqrt(13))*[3;2];
%F2 = (1/sqrt(13))*[3;2];

%CONV = 'c'; %convex energy density
%CONV = 'nc'; %non-convex energy density


%% Method to find minimum of E(v)
%solver = 'TWdirectFsolve';

%mu1=0;mu2=0;lambda=0;t1=0;t2=0; pv=0; regParam=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% p-Laplace problem

%% Problem Type
%problem = 'AbsValue_Square';
%problem = 'AbsValue_Lshape';
%problem = 'AbsValue_Square_exact';
%problem = 'AbsValue_Lshape_exact';

%probClass = 'AbsValue';

%% PDE Solver
%pdeSolver = 'AVP1'; 
%   DWRbool=0;

% Degree of Absolute Value: W(x)=|x|^p + regParam * |x|^2, Standard: |x|^2
%pv=2;
%regParam=0.0;

%% Method to find minimum of E(v)
%solver = 'AVdirectFsolve';

%F1=0;F2=0;CONV = '';
%mu1=0;mu2=0;lambda=0;t1=0;t2=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dual Calculation

if DWRbool
  %% Choose Goal Function
  % GOAL = 'GOAL_IM'; % integral middle: J(tau,v) = 1/|Omega1|  int_Omega1 v dx
                      % or               J(tau,v) = 1/|Omega1|  int_Omega1 Dv dx
  % GOAL = 'GOAL_DP'; % dual problem: J(tau,v) = int_Omega W*_epsilon(tau) dx
                      % or            J(tau,v) = int_Omega W*_epsilon(DW(Dv)) dx
  % GOAL = 'GOAL_PR'; % for rel-eff proof: J(tau,v) = ||DW(Dv)||^2

  %% Goal Functional 
  J = @functionalMeanIntegralAll;
  J4n = @functionalMeanIntegral4n;

  JuExact = SaveJuExact(problem,GOAL,pdeSolver);
  GetJuExact = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Maximum Number Degrees of Freedom
 maxNrDoF = 1000;

%% Mark Criterion
 mark = 'bulk';
% mark = 'uniform';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for the nonlinear solver 'fsolve' of MATLAB
%options = optimset('Display','iter','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',100,'TolFun',1e-10,'TolX',1e-10);
options = optimset('Display','off','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',10,'TolFun',1e-10,'TolX',1e-10);

%% Initialisation
p = initFFW(pdeSolver,problem,mark,maxNrDoF,probClass,solver,'redGreenBlue',estimate);

if DWRbool
    p.statics.J = J;
    p.statics.J4n = J4n;
    p.JuExact = JuExact; %or, if unknown, use Aitken at the end of this file for uniform refinement
    p.params.GOAL = GOAL;
end

p.params.options = options;
p.params.RHS = 'RHS2';
p.params.DWRbool = DWRbool;
p.problem.mu1 = mu1;
p.problem.mu2 = mu2;
p.problem.t1 = t1;
p.problem.t2 = t2;
p.problem.lambda = lambda;
p.problem.F1 = F1;
p.problem.F2 = F2;
p.params.CONV = CONV;
p.problem.pv=pv;
p.problem.regParam=regParam;
p.problem.probClass=probClass;

%% Compute Discrete Solution
%p = computeSolution(p);
p.params.rhsIntegtrateExactDegree = 19;
p.params.nonLinearExactIntegrateDegree = 5;
p = p.statics.run(p);

%% Show Discrete Solution
%figure(2)
%p.params.output.holdIt = false; 
%p.params.output.showIteratively = false;
%
%set(gcf,'Name','Displacement');
%p = show('drawU',p);
%view(15,20)



%% Show Grid
figure(nr)
p.params.output.holdIt = false; 

set(gcf,'Name','Grid');
p = show('drawGrid',p);

if DWRbool
  %% evaluation of the functional
  for curLvl = 2 : length(p.level)
    p = p.statics.J(p,curLvl);
  end

  if GetJuExact %get Aitken Extrapolation of J, you have to save the result in SaveJuExact.m
    p = aitkenExtrapolation('Ju',p);
    p.JuExact = p.aitkenExtrapolationJu;
    fprintf('\n Ju = %.15g\n',p.JuExact)
  end

  for curLvl = 2 : length(p.level)
    p.level(curLvl).JError = norm(p.level(curLvl).Ju-p.JuExact);
  end
end

if DrawSigma
 for curLvl = 2 : length(p.level)
%  fprintf('\n||sigma||^2=%.15g und ||sigma h||^2=%.15g',p.params.sigmaSqL2A,p.level(end).sigmaSqL2)
  p.level(curLvl).SigmaError = sqrt(abs(p.params.sigmaSqL2A-p.level(curLvl).sigmaSqL2));
 end
end

%figure(4);
%p = show('drawUexact',p);  %Discrete P1-solution bei ODRTGOAL

%% Show Convergence History
figure(1)
p.params.output.holdIt = true; 
p.params.output.showIteratively = false;
p.params.output.minDoF = 1;


if DWRbool
  p.params.output.name = [estimate, ' ', GOAL, ' ',mark, ' \eta'];
  p = show('drawError_estimatedError',p);
  %if loadField('p.level(end)','JError',p,false)
  %  p.params.output.name = [estimate, GOAL, mark, ' J(e)'];
  %  p = show('drawError_JError',p);
  %end
else
  p.params.output.name = [estimate, ' ', mark, ' \eta'];
  p = show('drawError_estimatedError',p);
end

if DrawSigma
  if BoolL43
    p.params.output.name = '||\sigma_h-\sigma||_{L^{4/3}(\Omega)}';
    p = show('drawError_L43errorPhminusP0',p);
    
%    p.params.square=1;
%    p.params.output.name = '||\sigma_h-\sigma||^2_{L^{4/3}(\Omega)}';
%    p = show('drawError_L43errorPhminusP0Square',p);
  else
    p.params.output.name = '||\sigma_h-\sigma||_{L^2(\Omega)}';
    p = show('drawError_L2errorPhminusP0',p);

%    p.params.square=1;
%    p.params.output.name = '||\sigma_h-\sigma||^2_{L^{2}(\Omega)}';
%    p = show('drawError_L2errorPhminusP0Square',p);
    
%    p.params.output.name = [estimate, mark, '||\sigma||^2_{L^2(\Omega)}-||\sigma_h||^2_{L^2(\Omega)}'];
%    p = show('drawError_SigmaError',p);
  end
end  
grid off


if DWRbool
  disp('Effectivity index');
  for curLvl = 1 : length(p.level)
    disp(p.level(curLvl).estimatedError/p.level(curLvl).JError);
  end
end

