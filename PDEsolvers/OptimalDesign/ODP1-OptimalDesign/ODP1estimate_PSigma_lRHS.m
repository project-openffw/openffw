function p = ODP1estimate_PSigma_lRHS(p)
%author: Lena Noack

%% INPUT
n4e = p.level(end).geom.n4e;
lvl = size(p.level,2);
e4ed = p.level(end).enum.e4ed;
curLvl = length(p.level);
length4ed = p.level(end).enum.length4ed;
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);
%% compute the error between p and sigma 

%for curEdge=1:length(e4ed) %||p_h - sigma_h||^2_L^2
%    if e4ed(curEdge,2) ~= 0
%       eta4Ed(curEdge) = integrate(n4e(e4ed(curEdge,1),:),lvl,degree,@integrand,p) + ...
%                         integrate(n4e(e4ed(curEdge,2),:),lvl,degree,@integrand,p);
%    else    
%       eta4Ed(curEdge) = integrate(n4e(e4ed(curEdge,1),:),lvl,degree,@integrand,p);
%    end
%end

eta4T = integrate(n4e,lvl,degree,@integrand,p); % ||p_h - sigma_h||^2_L^2(T)
%eta4ed = sqrt(eta4Ed); % ||p_h - sigma_h||_L^2

%% OUTPUT
p.level(end).etaT = sqrt(eta4T);                %||p_h - sigma_h||_L^2(T)
p.level(end).estimatedError = sum(sqrt(eta4T)); %||p_h - sigma_h||_L^2(Omega)
%p.level(end).etaEd = eta4ed;
%p.level(end).estimatedError = norm(eta4ed,2);

%% supply the integrand ||p_h - sigma_h||_L^2
function val = integrand(x,y,curElem,lvl,p)

grad_h = p.statics.grad_h;  % grad_h=p_h
sigma_h = p.statics.sigma_h; %sigma_h = DW(p_h)
midPoint4e = p.level(end).enum.midPoint4e;

curGradh = grad_h(midPoint4e(curElem,1),midPoint4e(curElem,1),curElem,lvl,p);
sigmah = sigma_h(x,y,curElem,lvl,p);

gradh = ones(length(x),1)*curGradh;

val = sum((gradh - sigmah).*(gradh - sigmah),2);
val = reshape(val,[1 1 length(x)]);
    