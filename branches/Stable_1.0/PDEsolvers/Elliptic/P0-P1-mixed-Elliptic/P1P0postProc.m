function p = P1P0postProc(p)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load geometry
n4e = p.level(end).geom.n4e;
c4n = p.level(end).geom.c4n;

% load discrete solution 
x = p.level(end).x;

% load enumerated data
nrNodes = p.level(end).nrNodes;
nrElems = p.level(end).nrElems;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = x(2*nrNodes+1:end)';

pU = x(1:nrNodes)';
pV = x(nrNodes+1:2*nrNodes)';

pU = pU(n4e);
pV = pV(n4e);

grad = zeros(3,2,nrElems);
grad(:,1,:) = pU';
grad(:,2,:) = pV';

%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.level(end).u = u;
p.level(end).grad = grad;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
