function midPoint4e = getMidPoint4e(c4n,n4e)

coord4n4e = reshape( c4n(n4e',:),3,[] );
s = sum(coord4n4e,1);
midPoint4e = reshape(s,[],2)/3;
