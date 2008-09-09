function val = matrixMultiplication(A,B)
%For given 3-dimensional matrices A ( dim(A) = [n m k] ) and 
%B ( dim(B) = [m l k] ) matrixMultiplication computes the elementwise
%matrix product A(k)*B(k)

dimMatrixA = size(A);
dimMatrixB = size(B);

if dimMatrixA(2) ~= dimMatrixB(1)
    error('The number of columns of matrix A must be equal with the number of rows of matrix B')
end

if ( length(dimMatrixA) > 2 && dimMatrixA(3) ~= dimMatrixB(3) )
    error('The third dimension must be equal for both matrices.')
end
    
permA = permute(A,[2 1 3]);
for k = 1:dimMatrixB(2)
    repA((k-1)*dimMatrixA(2)+1:k*dimMatrixA(2),:,:) = permA;
end
% repA = repmat(permA,[dimMatrixB(2) 1 1]);
linA = repA(:);

for k = 1:dimMatrixA(1)
    repB(:,(k-1)*dimMatrixB(2)+1:k*dimMatrixB(2),:) = B;
end
% repB = repmat(B,[1 dimMatrixA(1) 1]);
linB = repB(:);


val = linA.*linB;
val = reshape(val,dimMatrixA(2),[]);
val = sum(val,1)';
val = reshape(val,dimMatrixB(2),dimMatrixA(1),[]);
val = permute(val,[2 1 3]);
