function result = rowaddr(A,I,J)
% adresses A with the index matrices I and J point wise,
% i.e. result = A( I(j,k), J(j,k)).

% initialise result
result = zeros(size(I));

% if size of A is small enough we can use linear adressing
if(prod(size(A)) < 2^31)	
	for curCol = 1:size(I,2)
		linIndex = sub2ind(size(A),I(:,curCol),J(:,curCol));
		result(:,curCol) = full( A(linIndex) );
	end
else
	for j = 1:size(I,2)
		dummyI = I(:,j);
		dummyJ = J(:,j);
		dummyres = zeros(size(I,1),1);
		for k = 1:size(I,1)
			dummyres(k) = A(dummyI(k),dummyJ(k));
		end
		result(:,j) = dummyres;
	end
end
