function p = P1estimate(p)
%author: Lena Noack
%
%estimate.m estimate the energy error for a conforming P1-FE method 
%through computing the oscillations of f, and
%the jump of the flux. 
%
%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discrete gradient of solution
ph = p.level(end).grad4e;

% problem data
f = p.problem.f;

% load enumerated data
e4ed = p.level(end).enum.e4ed;
length4ed = p.level(end).enum.length4ed;
area4e = p.level(end).enum.area4e;
midpoint4e = p.level(end).enum.midPoint4e;
midPoint4ed = p.level(end).enum.midPoint4ed;
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
% Oscillations for edge patch
% boundary elements are considered to be their own neighbours
[I,dontUse] = find(e4ed == 0);
e4ed(I,2) = e4ed(I,1);

firstElem = e4ed(:,1);
secElem = e4ed(:,2);

area4firstElem = area4e(firstElem);
area4secElem = area4e(secElem);

f4firstElem = f(midpoint4e(firstElem,1),midpoint4e(firstElem,2),firstElem,end,p);
f4secElem = f(midpoint4e(secElem,1),midpoint4e(secElem,2),secElem,end,p);

f4patch = (area4firstElem.*f4firstElem + area4secElem.*f4secElem)./...
            (area4firstElem + area4secElem);

f4T1 = area4e(firstElem).*(f4firstElem - f4patch).^2;
f4T2 = area4e(secElem).*(f4secElem - f4patch).^2;

osc4ed = length4ed.^2 .* ( f4T1 + f4T2 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta4ed = sqrt(jump4ed + osc4ed);

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).jump4ed = sqrt(jump4ed);
p.level(end).osc4ed = sqrt(osc4ed);
p.level(end).etaEd = eta4ed;
p.level(end).estimatedError = norm(eta4ed,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%