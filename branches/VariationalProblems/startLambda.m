function p = start

%startDWR(estimator, window number for grid, draw (||sigma||^2 - ||sigma_h||^2)^0.5,lambda)
GOAL = 'GOAL_PR';
lambda = 0.0145;
estimate = 'estimateDeltaEnergy'; startDWR(estimate,3,1,GOAL,lambda)
estimate = 'estimate_AvgSigma'; startDWR(estimate,3,0,GOAL,lambda)
estimate = 'estimate_DWR'; startDWR(estimate,3,0,GOAL,lambda)

lambda = 0.0155;
estimate = 'estimateDeltaEnergy'; startDWR(estimate,3,1,GOAL,lambda)
estimate = 'estimate_AvgSigma'; startDWR(estimate,3,0,GOAL,lambda)
estimate = 'estimate_DWR'; startDWR(estimate,3,0,GOAL,lambda)

lambda = 0.016;
estimate = 'estimateDeltaEnergy'; startDWR(estimate,3,1,GOAL,lambda)
estimate = 'estimate_AvgSigma'; startDWR(estimate,3,0,GOAL,lambda)
estimate = 'estimate_DWR'; startDWR(estimate,3,0,GOAL,lambda)

lambda = 0.0165;
estimate = 'estimateDeltaEnergy'; startDWR(estimate,3,1,GOAL,lambda)
estimate = 'estimate_AvgSigma'; startDWR(estimate,3,0,GOAL,lambda)
estimate = 'estimate_DWR'; startDWR(estimate,3,0,GOAL,lambda)
