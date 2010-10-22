function p4n = getP4n(Db,DbEd,e4n,e4ed,n4e,nrNodes,nrElems)

% TODO initialize the vectors I and J

I = [];
J = [];

for k = 1:length(DbEd)
    nodes = Db(k,:);
    edge = DbEd(k);
    elem = e4ed(edge,1);
    freeNode = setdiff(n4e(elem,:),nodes);
    [dontUse,dontUse,neighbours] = find(e4n([nodes,freeNode],:));
    I = [I;freeNode*ones(length(unique(neighbours)),1)];
    J = [J;unique(neighbours)];
end

innerNodes = setdiff(1:nrNodes,[unique(I);unique(Db)]);

for k = 1:length(innerNodes)
    [dontUse,dontUse,neighbours] = find(e4n(innerNodes(k),:));
    I = [I;innerNodes(k)*ones(length(unique(neighbours)),1)];
    neighbours = neighbours(:);
    J = [J;unique(neighbours)];
end

p4n = sparse(I,J,ones(length(J),1),nrNodes,nrElems);