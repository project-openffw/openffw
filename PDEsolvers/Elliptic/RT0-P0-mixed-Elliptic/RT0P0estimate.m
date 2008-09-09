function p = RT0P0estimate(p)
%estimate.m estimate the energy error for the mixed RT0-P0-FE method 
%through computing the oscillations of f, and
%the jump of the flux. 
%
%authors: David Guenther, Jan Reininghaus

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load discrete solution
grad = p.level(end).grad;

% problem data
f = p.problem.f;

% load enumerated data
ed4e = p.level(end).enum.ed4e;
e4ed = p.level(end).enum.e4ed;
DbEd = p.level(end).enum.DbEd;
NbEd = p.level(end).enum.NbEd;
length4ed = p.level(end).enum.length4ed;
area4e = p.level(end).enum.area4e;
midpoint4e = p.level(end).enum.midPoint4e;
nrEdges = p.level(end).nrEdges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ed4e = ed4e';
I = ed4e(:);

phU = grad(:,1,:);
phV = grad(:,2,:);

S = phU(:);

[dontUse,index] = setdiff(I,[NbEd;DbEd]);
S(index) = -S(index);

U1 = accumarray(I,S,[nrEdges,1],@sum);

phUrot = phU([2 3 1],:);
S = phUrot(:);

S(index) = -S(index);
U2 = accumarray(I,S,[nrEdges,1],@sum);

S = phV(:);

S(index) = -S(index);
V1 = accumarray(I,S,[nrEdges,1],@sum);

phVrot = phV([2 3 1],:);
S = phVrot(:);

S(index) = -S(index);
V2 = accumarray(I,S,[nrEdges,1],@sum);

innerEdges = setdiff(1:nrEdges, [NbEd;DbEd]);
U1 = U1(innerEdges);
U2 = U2(innerEdges);
V1 = V1(innerEdges);
V2 = V2(innerEdges);

h = length4ed(innerEdges);
c1 = (1 - sqrt(1/3))/2;
c2 = (1 + sqrt(1/3))/2;

jump4ed = zeros(nrEdges,1);
jump4ed(innerEdges) = h.^2/2.*( (U1+c1*(U2-U1)).^2 + (U1+c2*(U2-U1)).^2 ...
							  + (V1+c1*(V2-V1)).^2 + (V1+c2*(V2-V1)).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oscillations 
% boundary elements are considered to be their own neighbours
[I,dontUse] = find(e4ed == 0);
e4ed(I,2) = e4ed(I,1);

firstElem = e4ed(:,1);
secElem = e4ed(:,2);

area4firstElem = area4e(firstElem);
area4secElem = area4e(secElem);

f4firstElem = f(midpoint4e(firstElem,1),midpoint4e(firstElem,2),p);
f4secElem = f(midpoint4e(secElem,1),midpoint4e(secElem,2),p);

f4patch = (area4firstElem.*f4firstElem + area4secElem.*f4secElem)./...
            (area4secElem + area4secElem);

f4T1 = area4e(firstElem).*(f4firstElem - f4patch).^2;
f4T2 = area4e(secElem).*(f4secElem - f4patch).^2;

osc4ed = ( length4ed.^2 .* ( f4T1 + f4T2 ) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta4ed = sqrt(jump4ed + osc4ed);

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).jump4ed = sqrt(jump4ed);
p.level(end).osc4ed = sqrt(osc4ed);
p.level(end).etaEd = eta4ed;
p.level(end).estimatedError = norm(eta4ed,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
