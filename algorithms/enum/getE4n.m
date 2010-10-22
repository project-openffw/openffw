function e4n = getE4n(n4e)

nrElems = size(n4e,1);
I = [n4e(:,1);n4e(:,2);n4e(:,3)];
J = [n4e(:,2);n4e(:,3);n4e(:,1)];
S = [1:nrElems,1:nrElems,1:nrElems];
e4n = sparse(I,J,S);
