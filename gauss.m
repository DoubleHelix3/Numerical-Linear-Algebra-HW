% Note: All of this code is very inefficient
%(it's also really terrible because I wrote this at 3am while half asleep)

A = [2 10 8 8 6; 1 4 -2 4 -1; 0 2 3 2 1; 3 8 3 10 9; 1 4 1 2 1];
[P,L,U] = gauss_partial_pivoting(A);

A_inv = gauss_inverse(P,L,U)

function A_inv = gauss_inverse(P,L,U)
    n = size(P,1);
    
    A_inv = zeros(n);
    I = eye(n);
    for k=1:n
        A_inv(:,k) = solve(P,L,U,I(:,k));
    end

end

%b = [52 14 12 51 15]';
%c = [50 4 12 48 12]';

%[P, L, U] = gauss_partial_pivoting(A);


%x_b = solve(P,L,U,b);
%x_c = solve(P,L,U,c);

% returns the index of the first maximal element in x
function i = max_index(x)
    n = size(x,1);
    i = 1;
    for k=1:n
        if abs(x(k)) > abs(x(i))
            i = k;
        end
    end
end

function A = swap_rows(A, i, j)
    temp = A(i,:);
    A(i,:) = A(j,:);
    A(j,:) = temp;
end

% Output is the PLU decomposition 
% PA = LU
function [P, L, U] = gauss_partial_pivoting(A)
    n = size(A,1);
    assert(n == size(A,1)); % ensure A is square

    P = eye(n);
    L = zeros(n);
    U = A;

    for k=1:(n-1)
        max = k - 1 + max_index(U(k:n,k));
        assert(A(max,max) ~= 0); % A should be non-singular

        U = swap_rows(U,max,k);
        P = swap_rows(P,max,k);
        L = swap_rows(L,max,k);

        for i=(k+1):n
            l = U(i,k)/U(k,k);
            L(i,k) = l;
            U(i,:) = U(i,:) - l*U(k,:);
        end
    end

    L = L + eye(n);
end

function x = solve(P,L,U,b)
    n = size(P,1);

    % permute the coordinates
    b = P * b;

    % Ly = b
    y = b;
    for k=1:n
        for i=1:(k-1)
            y(k) = y(k) - L(k,i)*y(i);
        end
        y(k) = y(k)/L(k,k);
    end

    % Ux = y
    x = y;
    for k=n:-1:1
        for i=k+1:n
            x(k) = x(k) - U(k,i)*x(i);
        end
        x(k) = x(k)/U(k,k);
    end
end
