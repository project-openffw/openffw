function p = ODP1estimate_PSigma_nlRHS(p)
%author: Lena Noack

%% INPUT
n4e = p.level(end).geom.n4e;
lvl = size(p.level,2);
e4ed = p.level(end).enum.e4ed;
curLvl = length(p.level);
degree = loadField('p.params','nonLinearExactIntegrateDegree',p,19);
%% compute the error between p and sigma 

%%%%% first part
grad_h = p.statics.grad_h;  % grad_h=p_h
e4ed = p.level(end).enum.e4ed;
length4ed = p.level(end).enum.length4ed;
midPoint4e = p.level(end).enum.midPoint4e;
nrElems = p.level(end).nrElems;
lvl = length(p.level);

normals4ed = p.level(lvl).enum.normals4ed;
nrEdges = p.level(end).nrEdges;

%% Compute the jump error of the stress
evalGrad = zeros(nrElems,2);
for curElem = 1:nrElems
    evalGrad(curElem,:) = grad_h(midPoint4e(curElem,1),midPoint4e(curElem,1),...
                                    curElem,lvl,p);
end

phU = evalGrad(:,1);
phV = evalGrad(:,2);

e4ed = sort(e4ed,2);
[I,dontUse] = find(e4ed == 0);
e4ed(I,1) = e4ed(I,2);

ph4edU1 = phU(e4ed(:,1));
ph4edU2 = phU(e4ed(:,2));
ph4edU = ph4edU1 - ph4edU2;

ph4edV1 = phV(e4ed(:,1));
ph4edV2 = phV(e4ed(:,2));
ph4edV = ph4edV1 - ph4edV2;

normal = zeros(nrEdges,2);
for curEdge = 1:nrEdges
    normal(curEdge,:) = normals4ed(curEdge,:);
end


%% ESTIMATE

jumpP4ed = length4ed.^2 .* ((ph4edU.*normal(:,1)).^2 + ...
           (ph4edV.*normal(:,2)).^2); %h_E^2 ([p_h]nu)^2
       
%%%% second part           
for curEdge=1:length(e4ed) %||p_h - sigma_h||^2_L^2
    if e4ed(curEdge,2) ~= 0
       eta4Ed(curEdge) = integrate(n4e(e4ed(curEdge,1),:),lvl,degree,@integrand,p) + ...
                         integrate(n4e(e4ed(curEdge,2),:),lvl,degree,@integrand,p);
    else    
       eta4Ed(curEdge) = integrate(n4e(e4ed(curEdge,1),:),lvl,degree,@integrand,p);
    end
end

eta4ed = sqrt(eta4Ed' + jumpP4ed); %(||p_h - sigma_h||^2_L^2 + h_E^2 ([p_h]nu)^2)^{1/2}

%% OUTPUT
p.level(end).etaEd = eta4ed;
p.level(end).estimatedError = norm(eta4ed,2);

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