format long

A = [8 1; -2 1];
q0 = [1; 1];

[eigenvectors,eigenvalues] = eig(A, "vector");
[lambda, index] = max(eigenvalues);
v = eigenvectors(:,index);

j = 0;
q = q0;
while norm(q - v) > 10^(-6)
    % iterate
    q_next = A*q;
    % rescaling
    s = norm(q_next);
    q_next = q_next/s;
    r = norm(q_next-v)/norm(q-v);

    j = j+1;
    q = q_next;
end

v
eigenvalues
