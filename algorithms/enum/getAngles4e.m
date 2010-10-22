function angles4e = getAngles4e(tangents4e)

tangent4eU = tangents4e(:,1,:);
tangent4eV = tangents4e(:,2,:);

dummyU = -tangent4eU .* tangent4eU([3 1 2],:,:);
dummyV = -tangent4eV .* tangent4eV([3 1 2],:,:);

dummy = dummyU + dummyV;
dummy = reshape(dummy,3,[]);
angles4e = acos(dummy)';
