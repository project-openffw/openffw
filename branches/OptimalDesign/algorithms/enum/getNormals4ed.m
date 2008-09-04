function normals4ed = getNormals4ed(normals4e,ed4e)

normals4e = permute(normals4e, [3 2 1 ] );
normals4ed(ed4e(:,1),:) = normals4e(:,:,1); 
normals4ed(ed4e(:,2),:) = normals4e(:,:,2);
normals4ed(ed4e(:,3),:) = normals4e(:,:,3);