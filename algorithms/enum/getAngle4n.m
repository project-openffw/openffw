function angle4n = getAngle4n(angles4e,n4e,nrElems,nrNodes)

angle4n = zeros(nrNodes,1);

for curElem = 1:nrElems
	curNodes = n4e(curElem,:);
	curAngles = angles4e(curElem,:);
	angle4n(curNodes) = angle4n(curNodes) + curAngles';
end
	
