function p = ODP1estimate_Jump(p)
%author: David Guenther, Lena Noack

%% INPUT 
sigma_h = p.statics.sigma_h;
e4ed = p.level(end).enum.e4ed;
length4ed = p.level(end).enum.length4ed;
midPoint4e = p.level(end).enum.midPoint4e;
nrElems = p.level(end).nrElems;
lvl = length(p.level);

%% Compute the jump error of the stress
evalSigma = zeros(nrElems,2);
for curElem = 1:nrElems
    evalSigma(curElem,:) = sigma_h(midPoint4e(curElem,1),midPoint4e(curElem,2),...
                                    curElem,lvl,p);
end

phU = evalSigma(:,1);
phV = evalSigma(:,2);

e4ed = sort(e4ed,2);
[I,dontUse] = find(e4ed == 0);
e4ed(I,1) = e4ed(I,2);

ph4edU1 = phU(e4ed(:,1));
ph4edU2 = phU(e4ed(:,2));
ph4edU = ph4edU1 - ph4edU2;

ph4edV1 = phV(e4ed(:,1));
ph4edV2 = phV(e4ed(:,2));
ph4edV = ph4edV1 - ph4edV2;

jump4ed = length4ed.^2 .* (ph4edU.^2 + ph4edV.^2);
eta4ed = jump4ed;

%% OUTPUT 
p.level(end).jump4ed = jump4ed;
p.level(end).etaEd = eta4ed;
p.level(end).estimatedError = norm(eta4ed,2).^2;
