function p = P3estimate(p)

% author: Joscha Gedicke

%% INPUT
length4ed = p.level(end).enum.length4ed;
n4e = p.level(end).geom.n4e;
n4ed = p.level(end).enum.n4ed;
ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = loadField('p.params','rhsIntegtrateExactDegree',p,3);

%% ESTIMATE
h_T = max(length4ed(ed4e)')';

nu_T = h_T.^2.*integrate(n4e,curLvl,degree,@funcHandleResiduum,p);
nu_E = length4ed.*integrate(n4ed,curLvl,3,@funcHandleNormalJump,p);

nu = sqrt( nu_T + 1/2*sum(nu_E(ed4e),2) );

%% OUTPUT
p.level(end).etaT = nu;
% p.level(end).etaOsc = nu;
% p.level(end).etaEd = 0;
p.level(end).estimatedError = norm(nu);