A = [4, 12, -16; 12, 37, -43; -16, -43, 98]
cholesky(A)

function A = cholesky(A)
    n = size(A,1);
    assert(size(A,2) == n) % check if square

    for i = 1:n
        for k = 1:i-1
            A(i,i) = A(i,i) - A(k,i)*A(k,i);
        end
        assert(A(i,i) > 0) % check if A positive definite
        A(i,i) = sqrt(A(i,i));
        
        for j = i+1:n
            for k = 1:i-1
                A(i,j) = A(i,j) - A(k,i)*A(k,j);
            end
            A(i,j) = A(i,j)/A(i,i);
        end
    end
end
