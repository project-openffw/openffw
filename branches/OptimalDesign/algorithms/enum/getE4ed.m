function e4ed = getE4ed(n4ed,e4n)

e4ed = rowaddr(e4n,n4ed,n4ed(:,[2 1]));
e4ed = sort(e4ed,2,'descend');