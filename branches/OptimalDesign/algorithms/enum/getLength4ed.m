function length4ed = getlength4ed(c4n,n4ed)

coordX = c4n(:,1);
coordY = c4n(:,2);

coordX4edge = coordX(n4ed);
coordY4edge = coordY(n4ed);

xLength = coordX4edge(:,2) - coordX4edge(:,1);
yLength = coordY4edge(:,2) - coordY4edge(:,1);

length4ed = sqrt(xLength.^2 + yLength.^2);
