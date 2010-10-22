function grad4e = getGrad4e(c4n,n4e,area4e)

xCoord4e = reshape(c4n(n4e,1),[],3);
yCoord4e = reshape(c4n(n4e,2),[],3);

grad4eY =  yCoord4e(:,[2,3,1]) - yCoord4e(:,[3,1,2]);
grad4eX =  xCoord4e(:,[3,1,2]) - xCoord4e(:,[2,3,1]);

grad4e = zeros(2*size(grad4eY,1),size(grad4eY,2));
grad4e(1:2:end,:) = grad4eY;
grad4e(2:2:end,:) = grad4eX;

dummy = repmat(area4e',2,1);
dummy = repmat(dummy(:),1,3);
grad4e = grad4e./(2*dummy);

grad4e = reshape(grad4e',3,2,[]);
