function p = ODP2estimate2(p)
% author: Lena Noack

%% INPUT
length4ed = p.level(end).enum.length4ed;
n4e = p.level(end).geom.n4e;
n4ed = p.level(end).enum.n4ed;
ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = loadField('p.params','rhsIntegrateExactDegree',p,5);

%% ESTIMATE
h_T = max(length4ed(ed4e)')';

etaT_A = h_T.^2.*integrate(n4e,curLvl,degree,@sigmaEquDiff,p);
eta_A = sqrt( etaT_A );

%% OUTPUT
p.level(end).etaT = eta_A;
p.level(end).estimatedError = norm(eta_A);



function val = sigmaEquDiff(x,y,curElem,lvl,p)

n4e = p.level(end).geom.n4e;
area4e = p.level(lvl).enum.area4e;
sigma_h = p.statics.sigma_h;

Sigma = sigma_h(x,y,curElem,lvl,p); 
n4e(curElem,:);
degree = loadField('p.params','rhsIntegrateExactDegree',p,5);
ASigma = (1/area4e(curElem))*integrate(n4e(curElem,:),lvl,degree,@SigmaV,p);


val(1,:,:) = norm(sum((Sigma - ASigma)'*(Sigma - ASigma),2));

function val = SigmaV(x,y,curElem,lvl,p)

sigma_h = p.statics.sigma_h;
val = norm((sigma_h(x,y,curElem,lvl,p))'*(sigma_h(x,y,curElem,lvl,p)));