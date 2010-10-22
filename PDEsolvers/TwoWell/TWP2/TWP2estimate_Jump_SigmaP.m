function p = TWP2estimate_Jump_SigmaP(p)
% Lena Noack


%% INPUT
length4ed = p.level(end).enum.length4ed;
n4e = p.level(end).geom.n4e;
n4ed = p.level(end).enum.n4ed;
ed4e = p.level(end).enum.ed4e;
curLvl = length(p.level);
degree = loadField('p.params','rhsIntegrateExactDegree',p,5);

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

jumpP4ed = length4ed.^(2) .* ((ph4edU.*normal(:,1)).^2 + ...
           (ph4edV.*normal(:,2)).^2); %h_E^2 ([p_h]nu)^2

jumpS4ed = length4ed.*sqrt(integrate(n4ed,curLvl,degree,@integrand,p)); %h_E ||[sigma_h]nu||^2_{L^4}

%eta4ed = sqrt(jumpS4ed).*sqrt(jumpP4ed); %h_E^{3/2} ||[sigma_h]nu||_{L^2} |[p_h]nu|
eta4ed = sqrt(sqrt(jumpS4ed).*sqrt(jumpP4ed)); %h_E^{3/2} ||[sigma_h]nu||_{L^4} |[p_h]nu|

%% OUTPUT 
p.level(end).jump4ed = sqrt(jumpS4ed); %h_E^0.5 ||[sigma_h]nu||_{L^2}
p.level(end).etaEd = eta4ed;
p.level(end).estimatedError = norm(eta4ed,2);


function val = integrand(x,y,curEdge,lvl,p)

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
val(1,1,:) = sum((evalSigma1 - evalSigma2).^4,2); %sum( (sigma_h|T_1 - sigma_h|T_2 ) * nu_E )^4
