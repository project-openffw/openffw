function p = ODP1prolong(p,lvl)
%author: David Guenther

%% INPUT
Db = p.level(lvl).geom.Db;
c4n = p.level(lvl).geom.c4n;
n4e = p.level(lvl).geom.n4e;
parents4e = p.level(lvl).enum.parents4e;
nrNodes = p.level(lvl).nrNodes;
nrElems = p.level(lvl).nrElems;

u_D = p.problem.u_D;
u_h = p.statics.u_h;

% x0 = p.level(1).P1x;

%% prolongation
DbNodes = unique(Db);
x0 = ones(nrNodes,1);

if lvl > 2
    x0 = zeros(nrNodes,1);
    for curElem = 1:nrElems
        nodes = n4e(curElem,:);
        coords = c4n(nodes,:);
        parent = parents4e(curElem);

        prolongU = u_h(coords(:,1),coords(:,2),parent,lvl-1,p);
        x0(nodes) = prolongU;
    end
end

x0(unique(Db)) = u_D(c4n(DbNodes,1),c4n(DbNodes,2),p);

%% OUTPUT
p.level(lvl).x = x0;
