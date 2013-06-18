function p = ODP1run(p)
%author: David Guenther

%% initialization
% Initial enumeration
enumerate = p.statics.enumerate;
p = enumerate(p);

% Set initial values
nrNodes = p.level(1).nrNodes;
p.level(1).x = zeros(nrNodes,1);

% supply the discrete functions
p.statics.basisU = @getDisplacementBasis;
p.statics.stressBasis = @getStressBasis;
p.statics.basis = @getBasis;

p = ODP1getUhFuncs(p);
p = ODP1getPhFuncs(p);

% compute the discrete solution
p = ODcomputeSolution(p);

%% supply displacement basis (P1)
function basisU = getDisplacementBasis(x,y,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
area4e = p.level(lvl).enum.area4e;

nodes = n4e(curElem,:);
coords = c4n(nodes,:);
area = area4e(curElem);

P1 = coords(1,:);
P2 = coords(2,:);
P3 = coords(3,:);

x = x';
y = y';

b1 = 1/2/area*( (P2(2)-P3(2))*x + (P3(1)-P2(1))*y + P2(1)*P3(2)-P3(1)*P2(2) );
b2 = 1/2/area*( (P3(2)-P1(2))*x + (P1(1)-P3(1))*y + P3(1)*P1(2)-P1(1)*P3(2) );
b3 = 1/2/area*( (P1(2)-P2(2))*x + (P2(1)-P1(1))*y + P1(1)*P2(2)-P2(1)*P1(2) );


basisU = [b1;b2;b3]';

%% supply stress basis
function basisSigma = getStressBasis(x,y,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
area4e = p.level(lvl).enum.area4e;

nodes = n4e(curElem,:);
coords = c4n(nodes,:);
area = area4e(curElem);

P1 = coords(1,:);
P2 = coords(2,:);
P3 = coords(3,:);

x = x';
y = y';

b1 = 1/2/area*[ (P2(2)-P3(2)), (P3(1)-P2(1)) ];
b2 = 1/2/area*[ (P3(2)-P1(2)), (P1(1)-P3(1)) ];
b3 = 1/2/area*[ (P1(2)-P2(2)), (P2(1)-P1(1)) ];

basisSigma = repmat([b1;b2;b3],[1,1,length(x)]);

function basisU = getBasis(x,y,curElem,lvl,p)

n4e = p.level(lvl).geom.n4e;
c4n = p.level(lvl).geom.c4n;
area4e = p.level(lvl).enum.area4e;

nodes = n4e(curElem,:);
coords = c4n(nodes,:);
area = area4e(curElem);

P1 = coords(1,:);
P2 = coords(2,:);
P3 = coords(3,:);

b1 = 1/2/area*( (P2(2)-P3(2))*x + (P3(1)-P2(1))*y + P2(1)*P3(2)-P3(1)*P2(2) );
b2 = 1/2/area*( (P3(2)-P1(2))*x + (P1(1)-P3(1))*y + P3(1)*P1(2)-P1(1)*P3(2) );
b3 = 1/2/area*( (P1(2)-P2(2))*x + (P2(1)-P1(1))*y + P1(1)*P2(2)-P2(1)*P1(2) );

basisU = zeros(size([x y],1),2,6);
basisU(:,1,1) = b1;
basisU(:,1,2) = b2;
basisU(:,1,3) = b3;
basisU(:,2,4) = b1;
basisU(:,2,5) = b2;
basisU(:,2,6) = b3;