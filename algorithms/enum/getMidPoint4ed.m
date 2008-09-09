function midpoint4ed = getMidPoint4ed(c4n,n4ed)

coordX = c4n(:,1);
coordY = c4n(:,2);

coordX4n4ed = coordX(n4ed);
coordY4n4ed = coordY(n4ed);

midpointX = sum(coordX4n4ed,2)/2;
midpointY = sum(coordY4n4ed,2)/2;

midpoint4ed = [midpointX,midpointY];
