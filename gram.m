m = 10;
V = hilb(m);

Q = gram_schmidt(V)
norm(eye(m) - Q'*Q)

% I know I could just get R from gram schmidt process
% but I'm too lazy
R = Q'*V;

norm(V - Q*R)

% gram schmidt
function Q = gram_schmidt(A)
    n = size(A,2);
    Q = zeros(n);
    for k=1:n
        Q(:,k) = A(:,k);

        for j=1:k-1
            Q(:,k) = Q(:,k) - (Q(:,j)' * A(:,k))*Q(:,j);
        end

        Q(:,k) = Q(:,k) / norm(Q(:,k));
    end
end
