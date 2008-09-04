function gradNC4e = getGradNC4e(midpoint4ed,ed4e,area4e)

xCoord4e = reshape(midpoint4ed(ed4e,1),[],3);
yCoord4e = reshape(midpoint4ed(ed4e,2),[],3);

gradNC4eY =  yCoord4e(:,[2,3,1]) - yCoord4e(:,[3,1,2]);
gradNC4eX =  xCoord4e(:,[3,1,2]) - xCoord4e(:,[2,3,1]);

gradNC4e = zeros(2*size(gradNC4eY,1),size(gradNC4eY,2));
gradNC4e(1:2:end,:) = gradNC4eY;
gradNC4e(2:2:end,:) = gradNC4eX;

dummy = repmat(area4e',2,1)/4;
dummy = repmat(dummy(:),1,3);
gradNC4e = gradNC4e./(2*dummy);

gradNC4e = reshape(gradNC4e',3,2,[]);
