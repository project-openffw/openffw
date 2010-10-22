function p = ODP1estimate_Jump_Avg(p)
%author: Lena Noack

%% INPUT
n4e = p.level(end).geom.n4e;
n4ed = p.level(end).enum.n4ed;
e4ed = p.level(end).enum.e4ed;
lvl = size(p.level,2);
length4ed = p.level(end).enum.length4ed;

degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);
%% compute the average error on patch T1 \cup T2 for T1 \cap T2 = E

for curEdge=1:length(e4ed) %||p_h-Ap_h||^2_{L^2}
    if e4ed(curEdge,2) ~= 0
       etaAvg(curEdge) = abs(integrate(n4e(e4ed(curEdge,1),:),lvl,degree,@integrand,p)) + ...
                         abs(integrate(n4e(e4ed(curEdge,2),:),lvl,degree,@integrand,p));
    else    
       etaAvg(curEdge) = abs(integrate(n4e(e4ed(curEdge,1),:),lvl,degree,@integrand,p));
    end
end

etaJump = length4ed.*integrate(n4ed,lvl,degree,@integrandJ,p); %h_E ||[sigma_h]nu||^2_{L^2}
eta4Ed = sqrt(sqrt(etaAvg)'.*sqrt(etaJump));

%% OUTPUT
p.level(end).etaEd = eta4Ed;
p.level(end).estimatedError = norm(eta4Ed,2);

%% supply the integrand ||p_h - Ap_h||_L^2
function val = integrand(x,y,curElem,lvl,p)

Ap_h = p.statics.Ap_h;
sigma_h = p.statics.sigma_h;

Aph = Ap_h(x,y,curElem,lvl,p);
sigmah = sigma_h(x,y,curElem,lvl,p);

val = sum((Aph - sigmah).*(Aph - sigmah),2);

val = reshape(val,[1 1 length(x)]);
    
function val = integrandJ(x,y,curEdge,lvl,p)

sigma_h = p.statics.sigma_h; %sigma_h = DW(p_h)
e4ed = p.level(lvl).enum.e4ed;
normals4ed = p.level(lvl).enum.normals4ed;

elems = e4ed(curEdge,:);
normal = normals4ed(curEdge,:);

evalSigma1 = sigma_h(x,y,elems(1),lvl,p)*normal'; %sigma_h|T_1 * nu_E

if elems(2) ~= 0
    % E is an interior edge
    evalSigma2 = sigma_h(x,y,elems(2),lvl,p)*normal'; %sigma_h|T_2 * nu_E
else
    % E is a Dirichlet Edge
    evalSigma2 = evalSigma1;
end

val = zeros(1,1,length(x));
val(1,1,:) = sum((evalSigma1 - evalSigma2).^2,2); %sum( (sigma_h|T_1 - sigma_h|T_2 ) * nu_E )^2
