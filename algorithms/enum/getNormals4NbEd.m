function normals4NbEd = getNormals4NbEd(c4n,Nb)

normals4NbEd = [];
if(~isempty(Nb))
	tangents = ( c4n(Nb(:,2),:) - c4n(Nb(:,1),:) );
	dummy = ones(size(tangents,1),2);
	normals = tangents(:,[2 1]).*[dummy(:,1),-dummy(:,2)];
	lengthNormal = sqrt(normals(:,1).^2 + normals(:,2).^2);

	normals4NbEd  = normals./[lengthNormal, lengthNormal];
end