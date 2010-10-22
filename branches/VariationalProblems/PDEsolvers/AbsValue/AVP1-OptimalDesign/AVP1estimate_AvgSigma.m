function p = P1estimate_Avg(p)
%author: David Guenther

%% INPUT
n4e = p.level(end).geom.n4e;
lvl = size(p.level,2);

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);
%% compute the average error 
eta4T = integrate(n4e,lvl,degree,@integrand,p);

%% OUTPUT
p.level(end).etaT = sqrt(eta4T);
p.level(end).estimatedError = sqrt(sum(eta4T));

%% supply the integrand ||p_h - Ap_h||_L^2
function val = integrand(x,y,curElem,lvl,p)

Ap_h = p.statics.Ap_h;
sigma_h = p.statics.sigma_h;

Aph = Ap_h(x,y,curElem,lvl,p);
sigmah = sigma_h(x,y,curElem,lvl,p);

val = sum((Aph - sigmah).*(Aph - sigmah),2);

val = reshape(val,[1 1 length(x)]);
    