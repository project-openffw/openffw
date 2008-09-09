function area4n = getArea4n(area4e,e4n)

[I,J,S] = find(e4n);
newS = area4e(S);
area4n = sparse(I,J,newS);
area4n = full(sum(area4n,2));
