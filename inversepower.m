A = [8 1; -2 1];
q0 = [1; 0];

[eigenvectors, eigenvalues] = eig(A, "vector")

[v,lambda] = rayliegh(A,q0, 4);

function [v,lambda] = rayliegh(A, q0, numiterations)
    [n,n] = size(A);

    q = q0;
    for j=1:numiterations
        prevq = q;
        p = q'*A*q/(q'*q); % rayliegh quotient
        q = (A - p*eye(n))\q;
        q = q/norm(q);
        
        q
    end

    % eigenvector
    v = q;
    % approximately the associated eigenvalue
    lambda = mean(A*v./v)
end
