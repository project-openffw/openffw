function ed4n = getEd4n(e4n)

ed4n = triu(e4n+e4n');
[I,J] = find(ed4n);
nrEdges = length(I);

ed4n = sparse(I,J,1:nrEdges,size(e4n,1),size(e4n,1));
ed4n = ed4n+ed4n';