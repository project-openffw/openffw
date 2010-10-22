function normals4DbEd = getNormals4DbEd(c4n,Db)

normals4DbEd = [];
if(~isempty(Db))
	tangents = ( c4n(Db(:,2),:) - c4n(Db(:,1),:) );
	dummy = ones(size(tangents,1),2);
	normals = tangents(:,[2 1]).*[dummy(:,1),-dummy(:,2)];
	lengthNormal = sqrt(normals(:,1).^2 + normals(:,2).^2);

	normals4DbEd  = normals./[lengthNormal, lengthNormal];
end
