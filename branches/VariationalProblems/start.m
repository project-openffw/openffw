function p = start

%startDWR(estimator, window number for grid, draw (||sigma||^2 - ||sigma_h||^2)^0.5)

estimate = 'estimateDeltaEnergy'; startDWR(estimate,3,0,'GOAL_DP')

%% eta - simple residual estimators
%estimate = 'estimate_Jump'; startDWR(estimate,4,0)
%estimate = 'estimate_Jump_SigmaP'; startDWR(estimate,5,0)
%estimate = 'estimate_AvgSigma'; startDWR(estimate,6,0)
%estimate = 'estimate_Avg'; startDWR(estimate,7,0)
%estimate = 'estimate_JumpCC'; startDWR(estimate,8,0)
%estimate = 'estimate_PD'; startDWR(estimate,9,0)
%estimate = 'estimate_Dual'; startDWR(estimate,10,0)

%% etaPWR - Estimators PWR
estimate = 'estimateH2'; startDWR(estimate,4,0)
estimate = 'estimateInterp'; startDWR(estimate,5,0)
estimate = 'estimateCompareSols'; startDWR(estimate,6,0)
estimate = 'estimate_PWR'; startDWR(estimate,7,0,'GOAL_IM')
estimate = 'estimate_PWR'; startDWR(estimate,8,0,'GOAL_DP')
estimate = 'estimate_PWR'; startDWR(estimate,9,0,'GOAL_PR')

%% etaDWR - Estimators DWR
%estimate = 'DWRestimateH2'; startDWR(estimate,4,0,'GOAL_PR')
%estimate = 'DWRestimateInterp'; startDWR(estimate,5,0,'GOAL_PR')
%estimate = 'DWRestimateCompareSols'; startDWR(estimate,6,0,'GOAL_PR')
%estimate = 'estimate_DWR'; startDWR(estimate,7,0,'GOAL_PR')









% estimateDeltaEnergy:      delta_ell = E(u_ell) - E(u)
% estimate_JumpSqrt:        h_E^{1/4} \|[\sigma_h]\nu\|^{1/2}_{L^2}
% estimate_Jump:            h_E^{1/2} \|[\sigma_h]\nu\|_{L^2}
% estimate_Jump_SigmaP:     h_E^{3/2} ||[sigma_h]nu||_{L^2} |[p_h]nu|
% estimate_AvgSigma:        ||sigma_h-A(sigma_h)||_L2
% estimate_Avg:             ||ph-A(ph)||_L2, ph=grad Uh
% estimate_JumpCC:          ( h_T^2 ||f_0 +div sigma_h||^2_(L^2)
%                             +0.5 sum h_E ||[sigma_h]nu||_{L^2} )^{1/2}
% estimate_PD:              E(u_ell) + int W*(tau) dx, div tau = -f
