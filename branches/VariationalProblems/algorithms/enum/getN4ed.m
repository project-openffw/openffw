function n4ed = getN4ed(ed4n)

% Create nodes4edge (n4ed) to quickly produce coordinates of the new nodes
% (ed4n is symetric)
[I,J,S] = find( triu(ed4n) );
n4ed( S,: ) = [I,J];
