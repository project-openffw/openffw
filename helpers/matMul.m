function val = matrixMultiplication(A,B)
%For given 3-dimensional matrices A ( dim(A) = [n m k] ) and
%B ( dim(B) = [m l k] ) matrixMultiplication computes the elementwise
%matrix product A(k)*B(k)

dimMatrixA = size(A);
dimMatrixB = size(B);

if ( length(dimMatrixA) > 2 && dimMatrixA(3) ~= dimMatrixB(3) )
    error('The third dimension must be equal for both matrices.')
end

if ( dimMatrixA(1) == 1 && dimMatrixA(2) == 1 )
    A1 = A(:)';
    A = zeros(dimMatrixB(1)*dimMatrixB(2),length(A1));
    for k = 1:dimMatrixB(1)*dimMatrixB(2)
      A(k,:) = A1;
    end
    % A = repmat(A,dimMatrixB(1)*dimMatrixB(2),1);
    A = reshape(A,[dimMatrixB(1),dimMatrixB(2),dimMatrixB(3)]);
    val = A.*B;
    return;
elseif ( dimMatrixB(1) == 1 && dimMatrixB(2) == 1 )
    B1 = B(:)';
    B = zeros(dimMatrixA(1)*dimMatrixA(2),length(B1));
    for k = 1:dimMatrixB(1)*dimMatrixB(2)
      B(k,:) = B1;
    end
    % B = repmat(B,dimMatrixA(1)*dimMatrixA(2),1);
    B = reshape(B,[dimMatrixA(1),dimMatrixA(2),dimMatrixA(3)]);
    val = B.*A;
    return;
else
    if dimMatrixA(2) ~= dimMatrixB(1)
        error('The number of columns of matrix A must be equal with the number of rows of matrix B')
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
end