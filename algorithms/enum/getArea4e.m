function area4e = getArea4e(c4n,n4e)

xCoord4e = reshape(c4n(n4e,1),[],3);
yCoord4e = reshape(c4n(n4e,2),[],3);

area4e = 1/2 * ( (xCoord4e(:,2)-xCoord4e(:,1)).*(yCoord4e(:,3)-yCoord4e(:,1))...
		- (yCoord4e(:,2)-yCoord4e(:,1)).*(xCoord4e(:,3)-xCoord4e(:,1)) );
