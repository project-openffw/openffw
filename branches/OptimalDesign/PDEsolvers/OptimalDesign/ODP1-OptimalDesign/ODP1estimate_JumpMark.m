function p = ODP1estimate_JumpMark(p)
%estimate.m estimate the energy error for a conforming P1-FE method 
%through computing the oscillations of f, and
%the jump of the flux. 
%
%authors: David Guenther, Jan Reininghaus

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discrete gradient of solution
ph = p.level(end).grad4e;

% problem data
f = p.problem.f;
% g = p.problem.g;

% load enumerated data
e4ed = p.level(end).enum.e4ed;
NbEd = p.level(end).enum.NbEd;
length4ed = p.level(end).enum.length4ed;
area4e = p.level(end).enum.area4e;
midpoint4e = p.level(end).enum.midPoint4e;
midPoint4ed = p.level(end).enum.midPoint4ed;
normals4NbEd = p.level(end).enum.normals4NbEd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phU = ph(:,1);
phV = ph(:,2);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta4ed = sqrt(jump4ed);


%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).jump4ed = sqrt(jump4ed);
p.level(end).etaEd = eta4ed;
p.level(end).estimatedError = norm(eta4ed,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


etaAvg = integrate(n4e,lvl,degree,@integrand,p);
p.level(end).etaAvg = sqrt(etaAvg);

function val = integrand(x,y,curElem,lvl,p)

Ap_h = p.statics.Ap_h;
sigma_h = p.statics.sigma_h;

Aph = Ap_h(x,y,curElem,lvl,p);
sigmah = sigma_h(x,y,curElem,lvl,p);

val = sum((Aph - sigmah).*(Aph - sigmah),2);

val = reshape(val,[1 1 length(x)]);
    