function p = startODoptimalLambda

%lambdaVals = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
lambdaVals = [-1,-0.1,0,0.1,0,1];
enVals = [];
for lambda = lambdaVals

pdeSolver = 'ODP1';
solver = 'ODdirectFsolveOL';
estimator = 'estimate_Jump';
problem = 'OptimalDesign_TwoSquare';
maxNrDoF = 1000;
minDoF = 10;
mark = 'bulk';
mu1 = 1;
mu2 = 2;

t1 = (2*lambda*mu1/mu2).^(1/2);
t2 = mu2/mu1*t1;

options = optimset('Display','off','Jacobian','on','NonlEqnAlgorithm','gn','MaxIter',10,'TolFun',1e-10,'TolX',1e-10);

%% COMPUTE DISCRETE SOLUTION
p = initFFW(pdeSolver,problem,mark,maxNrDoF,'optimalDesign',solver,'redGreenBlue',estimator);
p.params.options = options;
p.params.RHS = 'RHS2';
p.problem.mu1 = mu1;
p.problem.mu2 = mu2;

p.problem.t1 = t1;
p.problem.t2 = t2;
p.problem.lambda = lambda;

%dummy = (0.0 : 0.002 : 0.01).^3;
dummy = (0.01 : 0.125 : 1.01).^3;
p.params.epsSequence = dummy(end:-1:1);
p.params.epsPower = 1;
p.params.epsFactor = 1;
p.params.rhsIntegtrateExactDegree = 19;
p.params.nonLinearExactIntegrateDegree = 5;
p = p.statics.run(p);

energy = getEnergy(p);
enVals = [enVals,energy];


hold all
figure(4);
p.params.output.name = ['estError, \lambda=' num2str(lambda)];
p = show('drawError_estimatedError',p);
hold all
figure(5);
p.params.output.name = ['L2, \sigma-\sigma_h, \lambda=' num2str(lambda)];
p = show('drawError_L2errorPhminusP0',p);
hold all
figure(6);
p.params.output.name = ['L2, u-u_h, \lambda=' num2str(lambda)];
p = show('drawError_L2errorDisplacement',p);


end

%lambdaVals
%enVals

figure(7)
plot(lambdaVals',enVals')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% supply the discrete energy
function val = getEnergy(p)

energy_h = p.statics.energy_h;
lvl = size(p.level,2);
n4e = p.level(lvl).geom.n4e;

intVal = integrate(n4e,lvl,10,energy_h,p);
val = sum(intVal);
