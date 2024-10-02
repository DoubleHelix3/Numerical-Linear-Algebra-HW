% Generate A
num = 20; % or other appropriate number
G = numgrid('N',num);
B = delsq(G);
A = sprandn(B) + 1i*sprandn(B);
%spy(A)
[~,n] = size(A);

% inital vector
q = ones(n, 1);

m = 60; % = [15, 30, 45, 60];

H = arnoldi(A,q,m);
ritz = eig(H(1:m,1:m));

%spy(H)

hold off
eigval = eig(full(A));
plot (real(eigval),imag(eigval),'r+')
hold on
plot(real(ritz),imag(ritz),'bo')

function H = arnoldi(A, q, m)
    [~,n] = size(A);
    Q = zeros(n, m);
    Q(:,1) = q/norm(q);
    H = zeros(m+1, m);
    for k=1:m
        Q(:,k+1) = A*Q(:,k);
        H(1:k,k) = Q(:,1:k)'*Q(:,k+1);
        Q(:,k+1) = Q(:,k+1) - Q(:,1:k)*H(1:k,k); % orthogonalize
        s = Q(:,1:k)' * Q(:,k+1);
        Q(:,k+1) = Q(:,k+1) - Q(:, 1:k)*s; %reorthogonalize
        H(1:k,k) = H(1:k,k) + s;
        H(k+1,k) = norm(Q(:,k+1));

        %check to see if space is invarient
        tol = 10^(-6); % arbitrary
        assert(H(k+1,k) > tol)
    
        Q(:,k+1) = Q(:,k+1)/H(k+1,k);
    end
end
