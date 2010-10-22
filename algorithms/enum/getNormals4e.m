function normals4e = getNormals4e(c4n,n4e,length4ed,ed4e)

xCoord4e = reshape(c4n(n4e,1),[],3);
yCoord4e = reshape(c4n(n4e,2),[],3);

tangent4eX =  xCoord4e(:,[2,3,1]) - xCoord4e;
tangent4eY =  yCoord4e(:,[2,3,1]) - yCoord4e;

% Normalize 
lengthEd4e = reshape(length4ed(ed4e),[],3);
tangent4eX = tangent4eX ./ lengthEd4e;
tangent4eY = tangent4eY ./ lengthEd4e;

normals4e = zeros(2*size(tangent4eY,1),size(tangent4eY,2));
normals4e(1:2:end,:) = tangent4eY;
normals4e(2:2:end,:) = -tangent4eX;

normals4e = reshape(normals4e',3,2,[]);
