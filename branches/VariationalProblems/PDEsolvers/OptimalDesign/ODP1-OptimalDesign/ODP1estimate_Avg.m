function p = ODP1estimate_Avg(p)
%author: David Guenther

%% INPUT
n4e = p.level(end).geom.n4e;
lvl = size(p.level,2);

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);
%% compute the average error 
eta4T = integrate(n4e,lvl,degree,@integrand,p);

%% OUTPUT
p.level(end).etaT = sqrt(eta4T);
p.level(end).estimatedError = norm(eta4T,2);

%% supply the integrand ||p_h - Ap_h||_L^2
function val = integrand(x,y,curElem,lvl,p)

Ap_h = p.statics.AGrad_h;
grad_h = p.statics.grad_h;
%Ap_h = p.statics.Ap_h;
%sigma_h = p.statics.sigma_h;

Aph = Ap_h(x,y,curElem,lvl,p);
gradh = grad_h(x,y,curElem,lvl,p);

val = sum((Aph - gradh).*(Aph - gradh),2);

val = reshape(val,[1 1 length(x)]);
    